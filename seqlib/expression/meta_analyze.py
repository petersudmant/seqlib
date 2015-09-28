import argparse
from gtf_to_genes import *
import numpy as np
import scipy.stats as scp_stats
import pandas as pd

import logging
import pysam
import pdb
import math

class coverage_data():

    def __init__(self, g, typ = "longest_coding_transcript"):
        
        self.g = g
        self.strand = self.g.strand
        
        if typ == "longest_coding_transcript":
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

            self.exons = exons
            self.coding_exons = coding_exons
            
            self.l_UTR_s = exons[0][0] 
            self.l_UTR_e = coding_exons[0][0]
            self.l_UTR_l = self.l_UTR_e - self.l_UTR_s

            self.r_UTR_s = coding_exons[-1][1]
            self.r_UTR_e = exons[-1][1]
            self.r_UTR_l = self.r_UTR_e - self.r_UTR_s
            
            self.CDS_l = longest_t.coding_end - longest_t.coding_beg
            self.contig = "chr" in g.contig and g.contig or "chr%s"%g.contig
        else:
            assert True, "method %s not supported"%(typ)
        
        if g.strand:
            self.UTR_5p_l = self.l_UTR_l
            self.UTR_3p_l = self.r_UTR_l
            self.UTR_5p_s, self.UTR_5p_e = self.l_UTR_s, self.l_UTR_e 
            self.UTR_3p_s, self.UTR_3p_e = self.r_UTR_s, self.r_UTR_e 
        else:
            self.UTR_5p_l = self.r_UTR_l
            self.UTR_3p_l = self.l_UTR_l
            self.UTR_5p_s, self.UTR_5p_e = self.r_UTR_s, self.r_UTR_e 
            self.UTR_3p_s, self.UTR_3p_e = self.l_UTR_s, self.l_UTR_e 
        
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

    def get_seg_cvg(self, bamfile, contig, s, e):
        iter = bamfile.pileup(contig, s, e, truncate = True)
        cvg = np.zeros(e-s)
        for col in iter:
            n_non_del = np.sum(np.array([read.is_del==0 for read in col.pileups]))
            cvg[col.pos-s] = n_non_del
        return cvg
    
    def get_UTR_cvg(self, bamfile):
        self.UTR_5p_cvg = self.get_seg_cvg(bamfile, 
                                           self.contig, 
                                           self.UTR_5p_s, 
                                           self.UTR_5p_e)
        self.UTR_3p_cvg = self.get_seg_cvg(bamfile, 
                                           self.contig, 
                                           self.UTR_3p_s, 
                                           self.UTR_3p_e)

    def get_CDS_cvg(self, bamfile):
        cvgs = []
        for e in self.coding_exons:
            cvgs.append(self.get_seg_cvg(bamfile, 
                                         self.contig, 
                                         e[0], 
                                         e[1]))
        self.CDS_cvg = np.concatenate(cvgs)

    def get_cvg(self, bamfile):
        self.get_UTR_cvg(bamfile)
        self.get_CDS_cvg(bamfile)
        self.UTR_5p_mu  = np.mean(self.UTR_5p_cvg)
        self.UTR_3p_mu  =  np.mean(self.UTR_3p_cvg)
        self.CDS_mu = np.mean(self.CDS_cvg)
        self.total_mu = np.mean(np.concatenate([self.CDS_cvg, 
                                                self.UTR_5p_cvg, 
                                                self.UTR_3p_cvg]))
    
    def get_binned_cvg(self, vect, n_bins):
        
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
                "cds_len" : self.CDS_l,
                "UTR_5p_len" : self.UTR_5p_l,
                "UTR_3p_len" : self.UTR_3p_l,
                "cds_mu" : self.CDS_mu,
                "UTR_5p_mu" : self.UTR_5p_mu,
                "UTR_3p_mu" : self.UTR_3p_mu}
                

    def get_simple_summary_dict(self):
        dict = self.get_info_dict()
        
        dict.update({"UTR_5p_mu" : self.UTR_5p_mu,
                     "UTR_3p_mu" : self.UTR_3p_mu,
                     "CDS_mu" : self.UTR_3p_mu,
                     "mu_total" : self.total_mu})

        return dict

    def get_summary_dicts(self):
        UTR_5p = self.get_info_dict()
        UTR_3p = self.get_info_dict()
        CDS = self.get_info_dict()
        
        UTR_5p.update({"type" : "5p_UTR",
                       "pos" : 0,
                       "mu_depth" : self.UTR_5p_mu,
                       "mu_total" : self.total_mu})

        CDS.update({"type" : "CDS",
                    "pos" : 1,
                    "mu_depth" : self.CDS_mu,
                    "mu_total" : self.total_mu})
        
        UTR_3p.update({"type" : "3p_UTR",
                       "pos" : 2,
                       "mu_depth" : self.UTR_3p_mu,
                       "mu_total" : self.total_mu})
        
        return [UTR_5p, CDS, UTR_3p]                                                    
    

        
    def get_binned_cvg_dicts(self, n_bins):
        
        UTR_5p_cvg = self.get_binned_cvg(self.UTR_5p_cvg, n_bins)
        UTR_3p_cvg = self.get_binned_cvg(self.UTR_3p_cvg, n_bins)
        CDS_cvg = self.get_binned_cvg(self.CDS_cvg, n_bins)
        
        dicts = []
        
        for bin in range(n_bins): 
            UTR_5p = self.get_info_dict()
            UTR_3p = self.get_info_dict()
            CDS = self.get_info_dict()

            UTR_5p.update({"type" : "5p_UTR",
                           "bin" : bin, 
                           "pos" : 0, 
                           "cvg" : UTR_5p_cvg[bin],
                           "mu_total" : self.total_mu})

            CDS.update({"type" : "CDS",
                        "bin" : bin, 
                        "pos" : 1, 
                        "cvg" : CDS_cvg[bin],
                        "mu_total" : self.total_mu})

            UTR_3p.update({"type" : "3p_UTR",
                           "bin" : bin, 
                           "pos" : 2, 
                           "cvg" : UTR_3p_cvg[bin],
                           "mu_total" : self.total_mu})
            

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
        stop_codon_pos = stop_5p
        
        start_dist_sum = np.sum(START_DIST)
        stop_dist_sum = np.sum(STOP_DIST)

        for i in range(max([l_start_dist, l_stop_dist])):
            
            if i < l_start_dist:
                START_dict = self.get_info_dict()
                START_dict.update({"type" : "START_CODON",
                                   "pos" : i-start_codon_pos, 
                                   "cvg" : START_DIST[i],
                                   "sum_cvg" : start_dist_sum, 
                                   "strand" : self.strand, 
                                   "mu_total" : self.total_mu})
                dicts.append(START_dict)
            if i<l_stop_dist:
                STOP_dict = self.get_info_dict()
                STOP_dict.update({"type" : "STOP_CODON",
                                  "pos" : i-stop_codon_pos, 
                                  "cvg" : STOP_DIST[i],
                                  "sum_cvg" : stop_dist_sum, 
                                  "strand" : self.strand, 
                                  "mu_total" : self.total_mu})
                dicts.append(STOP_dict)
            
        return dicts


def get_cvg_objs_by_contig(contig_subset, genes, tr_contig, min_CDS, min_3p_UTR, min_5p_UTR):
    
    cvg_objs_by_contig = {}
    
    for g in genes['protein_coding']:
        contig = tr_contig(g.contig)
        if not contig in contig_subset:
            continue

        if not contig in cvg_objs_by_contig:
            cvg_objs_by_contig[contig] = []

        cvg_obj = coverage_data(g)
        if cvg_obj.pass_size_cutoff(min_CDS, min_3p_UTR, min_5p_UTR):
            cvg_objs_by_contig[contig].append(cvg_obj)

    for contig in cvg_objs_by_contig.keys():
        cvg_objs_by_contig[contig] = sorted(cvg_objs_by_contig[contig], key = lambda x: x.g.beg) 
    
    return cvg_objs_by_contig

    
if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_gtf_index", required=True)
    parser.add_argument("--fn_bam", required=True)
    parser.add_argument("--fn_out_summary", required=True)
    parser.add_argument("--fn_out_summary_simple", required=True)
    parser.add_argument("--fn_out_binned_cvg", required=True)
    parser.add_argument("--fn_out_CDS_start_stop", required=True)
    parser.add_argument("--gtf_ID", required=True)
    parser.add_argument("--contig_subset", default=None)
    parser.add_argument("--fn_logfile", default="/dev/null")
    parser.add_argument("--ENSEMBL_gtf", default=False, action="store_true")
    parser.add_argument("--n_cvg_bins", default=10, type=int)
    parser.add_argument("--min_CDS", default=1000, type=int)
    parser.add_argument("--min_3p_UTR", default=50, type=int)
    parser.add_argument("--min_5p_UTR", default=100, type=int)

    o = parser.parse_args()

    logger = logging.getLogger(o.fn_logfile)

    bamfile = pysam.AlignmentFile(o.fn_bam, 'rb')
    
    species_id, gtf_path, genes = get_indexed_genes_for_identifier(o.fn_gtf_index,
                                                                   logger, 
                                                                   o.gtf_ID)
    contig_subset = None
    if o.contig_subset: 
        contig_subset = o.contig_subset.split(":")
    else:
        contig_subset = ["chr%d"%i for i in range(1,24)]

    tr_contig = lambda x:  "chr" in x and x or "chr%s"%x
    
    cvg_objs_by_contig = get_cvg_objs_by_contig(contig_subset, 
                                                genes, 
                                                tr_contig,
                                                o.min_CDS,
                                                o.min_3p_UTR,
                                                o.min_5p_UTR)
    
    summary_outrows = []
    summary_simple_outrows = []
    binned_cvg_outrows = []
    CDS_start_stop_outrows = []

    for contig, cvg_objs in cvg_objs_by_contig.items():
        for i, cvg_obj in enumerate(cvg_objs):
            
            cvg_obj.get_cvg(bamfile)
            cvg_obj.print_summary()
            
            summary_outrows.extend(cvg_obj.get_summary_dicts())
            summary_simple_outrows.extend(cvg_obj.get_simple_summary_dict())
            binned_cvg_outrows.extend(cvg_obj.get_binned_cvg_dicts(o.n_cvg_bins))
            CDS_start_stop_outrows.extend(cvg_obj.get_CDS_start_stop_dicts())
            if i>1000:
                break
        break

    T_summary = pd.DataFrame(summary_outrows)
    T_summary.to_csv(o.fn_out_summary, index=False, sep="\t")
    
    T_summary_simple = pd.DataFrame(summary_simple_outrows)
    T_summary.to_csv(o.fn_out_summary_simple, index=False, sep="\t")

    T_binned_cvg = pd.DataFrame(binned_cvg_outrows)
    T_binned_cvg.to_csv(o.fn_out_binned_cvg, index=False, sep="\t")
    
    T_CDS_start_stop = pd.DataFrame(CDS_start_stop_outrows)
    T_CDS_start_stop.to_csv(o.fn_out_CDS_start_stop, index=False, sep="\t")





