import argparse
from gtf_to_genes import *
import numpy as np
import pandas as pd

import logging
import pysam
import pdb
import math

class coverage_data():

    def __init__(self, g, typ = "longest_coding_transcript"):
        
        self.g = g
        
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
        iter = bamfile.pileup(contig, s, e)
        cvg = np.array([col.nsegments for col in iter])
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
    
    def print_summary(self):
        pattern = ("Gene: {gene_name} {contig}:{start}-{end}\n"
                   "\t5' UTR mean: {utr5_mu}\n"
                   "\t3' UTR mean: {utr3_mu}\n"
                   "\tCDS mean: {cds_mu}\n")

        pattern = pattern.format(gene_name = self.g.names[0],
                                 contig = self.contig,
                                 start = self.g.beg,
                                 end = self.g.end,
                                 utr5_mu = self.UTR_5p_mu,
                                 utr3_mu = self.UTR_3p_mu,
                                 cds_mu = self.CDS_mu)
        print(pattern)
    
    def get_info_dict(self):
        return {"gene" : self.g.names[0],
                "contig" : self.contig,
                "start" : self.g.beg,
                "end" : self.g.end,
                "cds_len" : self.CDS_l,
                "UTR_5p_len" : self.UTR_5p_l,
                "UTR_3p_len" : self.UTR_3p_l}

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
        
        for bin in range(n_bins): 
            UTR_5p = self.get_info_dict()
            UTR_3p = self.get_info_dict()
            CDS = self.get_info_dict()
             
        return [{"gene" : self.g.names[0],
                "contig" : self.contig,
                "start" : self.g.beg,
                "end" : self.g.end,
                "cds_len" : self.CDS_l,
                "UTR_5p_len" : self.UTR_5p_l,
                "UTR_3p_len" : self.UTR_3p_l,
                "UTR_5p_mu" : np.mean(self.UTR_5p_cvg),
                "UTR_3p_mu" : np.mean(self.UTR_3p_cvg),
                "CDS_mu" : np.mean(self.CDS_cvg)}]


    def get_CDS_start_stop_dicts(self):

        return [{"gene" : self.g.names[0],
                "contig" : self.contig,
                "start" : self.g.beg,
                "end" : self.g.end,
                "cds_len" : self.CDS_l,
                "UTR_5p_len" : self.UTR_5p_l,
                "UTR_3p_len" : self.UTR_3p_l,
                "UTR_5p_mu" : np.mean(self.UTR_5p_cvg),
                "UTR_3p_mu" : np.mean(self.UTR_3p_cvg),
                "CDS_mu" : np.mean(self.CDS_cvg)}]

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
    parser.add_argument("--fn_out_binned_cvg", required=True)
    parser.add_argument("--fn_out_CDS_start_stop", required=True)
    parser.add_argument("--gtf_ID", required=True)
    parser.add_argument("--contig_subset", default=None)
    parser.add_argument("--fn_logfile", default="/dev/null")
    parser.add_argument("--ENSEMBL_gtf", default=False, action="store_true")
    parser.add_argument("--cvg_bins", default=1000, type=int)
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
    binned_cvg_outrows = []
    CDS_start_stop_outrows = []

    for contig, cvg_objs in cvg_objs_by_contig.items():
        for i, cvg_obj in enumerate(cvg_objs):
            cvg_obj.get_cvg(bamfile)
            cvg_obj.print_summary()
            
            summary_outrows.extend(cvg_obj.get_summary_dicts())
            binned_cvg_outrows.extend(cvg_obj.get_binned_cvg_dicts())
            CDS_start_stop_outrows.extend(cvg_obj.get_CDS_start_stop_dicts())
            if i>20:
                break
        break

    T_summary = pd.DataFrame(summary_outrows)
    T_summary.to_csv(o.fn_out_summary, index=False, sep="\t")

    T_binned_cvg = pd.DataFrame(binned_cvg_outrows)
    T_binned_cvg.to_csv(o.fn_out_binned_cvg, index=False, sep="\t")
    
    T_CDS_start_stop = pd.DataFrame(CDS_start_stop_outrows)
    T_CDS_start_stop.to_csv(o.fn_out_CDS_start_stop, index=False, sep="\t")





