import numpy as np
import scipy.stats as scp_stats
import pandas as pd

import pysam_ext.pysam_ext as pysam_ext
import itertools
import time
import pdb


def longest_coding_t(g):

    longest_t = None
    longest_l = None
    for t in g.transcripts:
        if longest_t == None and t.coding_beg != None:
            longest_t = t
            longest_t_l = np.absolute(t.coding_end - t.coding_beg)
        elif t.coding_end!=None and \
                 longest_t_l!=None and \
                 np.absolute(t.coding_end - t.coding_beg) > longest_t_l: 
            longest_t = t
            longest_t_l = np.absolute(longest_t.coding_end - longest_t.coding_beg) 
    
    exons = sorted(longest_t.get_exons())
    coding_exons = sorted(longest_t.get_coding_exons())

    contig = "chr" in g.contig and g.contig or "chr%s"%g.contig
    
    return exons, coding_exons, contig

def get_l_r_UTR_exons(exons, coding_exons):
    """
    NEED TO MAKE l_UTR_exons and r_UTR_exons
    """
    coding_ex_l, coding_ex_r = coding_exons[0], coding_exons[-1]
    
    l_UTR_exons = []
    r_UTR_exons = []
    
    for ex in exons:
        if ex[1]>coding_ex_l[0]:
            l_UTR_exons.append([ex[0],coding_ex_l[0]])
            break
        else:
            l_UTR_exons.append(ex)
    
    for ex in exons[::-1]:
        if ex[0]<coding_ex_r[1]:
            r_UTR_exons.append([coding_ex_r[1],ex[1]])
            break
        else:
            r_UTR_exons.append(ex)

    return l_UTR_exons, r_UTR_exons

class CoverageData():

    def __init__(self, g, typ = "virtual_gene"):
        
        self.g = g
        self.strand = self.g.strand
        
        if typ == "virtual_gene":
            self.contig = "chr" in g.contig and g.contig or "chr%s"%g.contig
            self.exons = sorted(g.virtual_exons)
            self.coding_exons = sorted(g.virtual_coding_exons)
        elif typ == "longest_coding_transcript":
            self.exons, self.coding_exons, self.contig = longest_coding_t(g)
        else:
            assert True, "method %s not supported"%(typ)
        
        l_UTR_exons, r_UTR_exons = get_l_r_UTR_exons(self.exons, self.coding_exons)

        if g.strand:
            self.UTR_5p_exons = sorted(l_UTR_exons)
            self.UTR_3p_exons = sorted(r_UTR_exons)
        else:
            self.UTR_5p_exons = sorted(r_UTR_exons)
            self.UTR_3p_exons = sorted(l_UTR_exons)
        
        self.CDS_l = np.sum(np.array([ex[1]-ex[0] for ex in self.coding_exons]))
        self.UTR_5p_l = np.sum(np.array([ex[1]-ex[0] for ex in self.UTR_5p_exons]))
        self.UTR_3p_l = np.sum(np.array([ex[1]-ex[0] for ex in self.UTR_3p_exons]))
        
        #RNA - space
        self.RNAcoord_UTR_5p_exons = []
        self.RNAcoord_UTR_3p_exons = []
        self.RNAcoord_coding_exons = []
        for i,e in enumerate(self.UTR_5p_exons):
            elen = e[1]-e[0]
            last_e = i>0 and self.RNAcoord_UTR_5p_exons[-1][1] or 0
            self.RNAcoord_UTR_5p_exons.append((last_e, last_e+elen))
            
        for i,e in enumerate(self.UTR_3p_exons):
            elen = e[1]-e[0]
            last_e = i>0 and self.RNAcoord_UTR_3p_exons[-1][1] or 0
            self.RNAcoord_UTR_3p_exons.append((last_e, last_e+elen))
        
        for i,e in enumerate(self.coding_exons):
            elen = e[1]-e[0]
            last_e = i>0 and self.RNAcoord_coding_exons[-1][1] or 0
            self.RNAcoord_coding_exons.append((last_e, last_e+elen))
        
        self.UTR_5p_cvg = None
        self.UTR_3p_cvg = None
        self.CDS_cvg = None

    def pass_size_cutoff(self, min_CDS, min_3p, min_5p): 
        if self.CDS_l >= min_CDS and \
                self.UTR_5p_l >= min_5p and \
                self.UTR_3p_l >= min_3p:
            return True
        else:
            return False

    def c_get_seg_cvg(self, bamfile, contig, s, e):
        iter = bamfile.pileup(contig, s, e, truncate = True)
        cvg = pysam_ext.get_seg_cvg(iter, e-s)
        return cvg

    def get_seg_cvg(self, bamfile, contig, s, e):
        return self.c_get_seg_cvg(bamfile, contig, s, e)
        """
        iter = bamfile.pileup(contig, s, e, truncate = True)
        cvg = np.zeros(e-s, dtype=np.int)
        for col in iter:
            n_non_del = np.sum(np.array([read.is_del==0 for read in col.pileups]))
            cvg[col.pos-s] = n_non_del
        print("py: %d %f"%(np.sum(cvg),t))
        return cvg
        """
    
    def get_cvg_over_loci(self, bamfile, loci):
        cvgs = []
        for loc in loci:
            cvgs.append(self.get_seg_cvg(bamfile, 
                                         self.contig, 
                                         loc[0], 
                                         loc[1]))
        return np.concatenate(cvgs)
    

    def get_cvg(self, bamfile):
        
        self.UTR_5p_cvg = self.get_cvg_over_loci(bamfile, self.UTR_5p_exons)
        self.UTR_3p_cvg = self.get_cvg_over_loci(bamfile, self.UTR_3p_exons)
        self.CDS_cvg = self.get_cvg_over_loci(bamfile, self.coding_exons)
        
        self.UTR_5p_mu  = np.mean(self.UTR_5p_cvg)
        self.UTR_3p_mu  =  np.mean(self.UTR_3p_cvg)
        self.CDS_mu = np.mean(self.CDS_cvg)
        
        self.UTR_5p_median  = np.median(self.UTR_5p_cvg)
        self.UTR_3p_median  =  np.median(self.UTR_3p_cvg)
        self.CDS_median = np.median(self.CDS_cvg)
        
        cat_cvg = np.concatenate([self.CDS_cvg, 
                                  self.UTR_5p_cvg, 
                                  self.UTR_3p_cvg])

        self.GENE_mu = np.mean(cat_cvg)
        self.GENE_median = np.median(cat_cvg)
        
    def get_binned_cvg(self, vect, n_bins):
        """
        Returning the mean
        """
        return scp_stats.binned_statistic(np.arange(vect.shape[0]),
                                          vect, 
                                          bins = n_bins)[0]

    def print_summary(self):
        pattern = ("Gene: {gene_name} {contig}:{start}-{end}\n"
                   "\t5' ({utr5_len}bp) UTR mean: {utr5_mu}\n"
                   "\t3' ({utr3_len}bp) UTR mean: {utr3_mu}\n"
                   "\tCDS ({cds_len}bp) mean: {cds_mu}\n")

        pattern = pattern.format(gene_name = self.g.names[0],
                                 contig = self.contig,
                                 start = self.g.beg,
                                 end = self.g.end,
                                 utr5_mu = self.UTR_5p_mu,
                                 utr3_mu = self.UTR_3p_mu,
                                 cds_mu = self.CDS_mu,
                                 utr5_len = self.UTR_5p_l,
                                 utr3_len = self.UTR_3p_l,
                                 cds_len = self.CDS_l)
        print(pattern)
    
    def get_info_dict(self):
        return {"gene" : self.g.names[0],
                "contig" : self.contig,
                "start" : self.g.beg,
                "end" : self.g.end,
                "CDS_len" : self.CDS_l,
                "UTR_5p_len" : self.UTR_5p_l,
                "UTR_3p_len" : self.UTR_3p_l,
                "CDS_mu" : self.CDS_mu,
                "UTR_5p_mu" : self.UTR_5p_mu,
                "UTR_3p_mu" : self.UTR_3p_mu,
                "CDS_median" : self.CDS_median,
                "UTR_5p_median" : self.UTR_5p_median,
                "UTR_3p_median" : self.UTR_3p_median, 
                "GENE_mu":self.GENE_mu,
                "GENE_median":self.GENE_median}
                

    def get_simple_summary_dict(self):
        dict = self.get_info_dict()
        return dict
    
    def get_bp_cvg_dicts(self):
            
        dicts = []
        
        UTR_5p_poses = list(itertools.chain(*[range(e[0],e[1]) for e in self.UTR_5p_exons]))
        UTR_3p_poses = list(itertools.chain(*[range(e[0],e[1]) for e in self.UTR_3p_exons]))
        cds_exon_poses = list(itertools.chain(*[range(e[0],e[1]) for e in self.coding_exons]))
        
        for i in xrange(self.UTR_5p_cvg.shape[0]):
            if self.UTR_5p_cvg[i]!=0:
                d = self.get_info_dict()
                d.update({"type" : "5p_UTR",
                           "pos" : UTR_5p_poses[i],
                           "cvg" : self.UTR_5p_cvg[i]})
                dicts.append(d)

        for i in xrange(self.UTR_3p_cvg.shape[0]):
            if self.UTR_3p_cvg[i]!=0:
                d = self.get_info_dict()
                d.update({"type" : "3p_UTR",
                           "pos" : UTR_3p_poses[i],
                           "cvg" : self.UTR_3p_cvg[i]})
                dicts.append(d)

        for i in xrange(self.CDS_cvg.shape[0]):
            if self.CDS_cvg[i]!=0:
                d = self.get_info_dict()
                d.update({"type" : "CDS",
                           "pos" : cds_exon_poses[i],
                           "cvg" : self.CDS_cvg[i]})
                dicts.append(d)
        
        """
        PUT 0s on the edges of exons - NOTE: not showing the intron loci
        even though they could have cvg
        """
        for e in self.UTR_5p_exons+self.UTR_3p_exons+self.coding_exons:
            d1 = self.get_info_dict()
            d2 = self.get_info_dict()
            
            d1.update({"type" : "edge",
                       "pos" : e[0]-1,
                       "cvg" : 0})

            d2.update({"type" : "edge",
                       "pos" : e[1]+1,
                       "cvg" : 0})
            dicts.append(d1)
            dicts.append(d2)
        return dicts

    def get_by_exon_dicts(self):
        return_dicts = []
        for i, UTR_5p_e in enumerate(self.UTR_5p_exons):
            d = self.get_info_dict()
            s,e  = UTR_5p_e
            s_i, e_i = self.RNAcoord_UTR_5p_exons[i]
            
            mu =  np.mean(self.UTR_5p_cvg[s_i:e_i])
            med = np.median(self.UTR_5p_cvg[s_i:e_i])

            d.update({"type": "5p_UTR",
                      "exon":i,
                      "exon_start":s,
                      "exon_end":e,
                      "mu_cvg": mu, 
                      "median_cvg": med})
            return_dicts.append(d)
        
        for i, CDS_e in enumerate(self.coding_exons):
            d = self.get_info_dict()
            s,e  = CDS_e
            s_i, e_i = self.RNAcoord_coding_exons[i]

            mu =  np.mean(self.CDS_cvg[s_i:e_i])
            med = np.median(self.CDS_cvg[s_i:e_i])

            d.update({"type": "CDS",
                      "exon":i,
                      "exon_start":s,
                      "exon_end":e,
                      "mu_cvg": mu,
                      "median_cvg": med})
            return_dicts.append(d)
        
        for i, UTR_3p_e in enumerate(self.UTR_3p_exons):
            d = self.get_info_dict()
            s,e  = UTR_3p_e
            s_i, e_i = self.RNAcoord_UTR_3p_exons[i]
            
            mu =  np.mean(self.UTR_3p_cvg[s_i:e_i])
            med = np.median(self.UTR_3p_cvg[s_i:e_i])
            
            d.update({"type": "3p_UTR",
                      "exon":i,
                      "exon_start":s,
                      "exon_end":e,
                      "mu_cvg": mu,
                      "median_cvg": med})
            return_dicts.append(d)

        return return_dicts

    def get_summary_dicts(self):
        UTR_5p = self.get_info_dict()
        UTR_3p = self.get_info_dict()
        CDS = self.get_info_dict()
        
        UTR_5p.update({"type" : "5p_UTR",
                       "pos" : 0,
                       "mu_cvg" : self.UTR_5p_mu})

        CDS.update({"type" : "CDS",
                    "pos" : 1,
                    "mu_cvg" : self.CDS_mu})
        
        UTR_3p.update({"type" : "3p_UTR",
                       "pos" : 2,
                       "mu_cvg" : self.UTR_3p_mu})
        
        return [UTR_5p, CDS, UTR_3p]                                                    
    
    def get_binned_cvg_dicts(self, n_bins):
        
        b_UTR_5p_cvg = self.get_binned_cvg(self.UTR_5p_cvg, n_bins)
        b_UTR_3p_cvg = self.get_binned_cvg(self.UTR_3p_cvg, n_bins)
        b_CDS_cvg = self.get_binned_cvg(self.CDS_cvg, n_bins)
        
        dicts = []
        
        for bin in range(n_bins): 
            UTR_5p = self.get_info_dict()
            UTR_3p = self.get_info_dict()
            CDS = self.get_info_dict()

            UTR_5p.update({"type" : "5p_UTR",
                           "bin" : bin, 
                           "pos" : 0, 
                           "cvg" : b_UTR_5p_cvg[bin]})

            CDS.update({"type" : "CDS",
                        "bin" : bin, 
                        "pos" : 1, 
                        "cvg" : b_CDS_cvg[bin]})

            UTR_3p.update({"type" : "3p_UTR",
                           "bin" : bin, 
                           "pos" : 2, 
                           "cvg" : b_UTR_3p_cvg[bin]})

            dicts.extend([UTR_5p, CDS, UTR_3p])
            
        return dicts


    def get_CDS_start_stop_dicts(self):
        """
        START CODON -50 3 +100
        STOP CODON -100 3 +100
        """
        start_3p = 50
        start_5p = 100
        
        stop_3p = 100
        stop_5p = 100

        START_DIST = None
        STOP_DIST = None
        
        if self.strand:
            """FWD"""
            START_DIST = np.concatenate([self.UTR_5p_cvg[-start_3p::],
                                         self.CDS_cvg[:start_5p+3]])
            STOP_DIST = np.concatenate([self.CDS_cvg[-stop_3p-3::],
                                        self.UTR_3p_cvg[:stop_5p]])
        else:
            """REV"""
            START_DIST = np.concatenate([self.CDS_cvg[-start_5p-3::],
                                         self.UTR_5p_cvg[:start_3p]])
                                         
            STOP_DIST = np.concatenate([self.UTR_3p_cvg[-stop_5p::],
                                        self.CDS_cvg[:stop_3p+3]])
            START_DIST = START_DIST[::-1]
            STOP_DIST = STOP_DIST[::-1]
        
        dicts = []
        l_start_dist = START_DIST.shape[0]
        l_stop_dist = STOP_DIST.shape[0]
        start_codon_pos = start_3p
        stop_codon_pos = stop_3p
        
        start_dist_sum = np.sum(START_DIST)
        stop_dist_sum = np.sum(STOP_DIST)

        for i in range(max([l_start_dist, l_stop_dist])):
            
            if i < l_start_dist:
                START_dict = self.get_info_dict()
                START_dict.update({"type" : "START_CODON",
                                   "pos" : i-start_codon_pos, 
                                   "cvg" : START_DIST[i],
                                   "sum_cvg" : start_dist_sum, 
                                   "strand" : self.strand })
                dicts.append(START_dict)
            if i<l_stop_dist:
                STOP_dict = self.get_info_dict()
                STOP_dict.update({"type" : "STOP_CODON",
                                  "pos" : i-stop_codon_pos, 
                                  "cvg" : STOP_DIST[i],
                                  "sum_cvg" : stop_dist_sum, 
                                  "strand" : self.strand })
                dicts.append(STOP_dict)
            
        return dicts
