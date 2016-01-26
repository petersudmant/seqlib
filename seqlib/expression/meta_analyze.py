import argparse
from gtf_to_genes import *
import numpy as np
import scipy.stats as scp_stats
import pandas as pd

import pysam_ext.pysam_ext as pysam_ext
from coveragedata import CoverageData 

import logging
import pysam
import pdb
import math
import time


def get_cvg_objs_by_contig(contig_subset, indiv_gene, genes, tr_contig, min_CDS, min_3p_UTR, min_5p_UTR):
    
    cvg_objs_by_contig = {}
    for g in genes['protein_coding']:
        contig = tr_contig(g.contig)
        if not contig in contig_subset:
            continue
        if indiv_gene and not indiv_gene in g.names:
            continue
        
        if not contig in cvg_objs_by_contig:
            cvg_objs_by_contig[contig] = []

        cvg_obj = CoverageData(g)
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
    parser.add_argument("--fn_out_summary_by_exon", required=True)

    parser.add_argument("--fn_out_binned_cvg", required=False)
    parser.add_argument("--fn_out_CDS_start_stop", required=False)
    
    parser.add_argument("--CDS_start_stop", default=False, action="store_true")
    parser.add_argument("--binned_cvg", default=False, action="store_true")
    
    parser.add_argument("--individual_gene", required=False, default=None)
    parser.add_argument("--gtf_ID", required=True)
    parser.add_argument("--contig_subset", default=None)
    parser.add_argument("--fn_logfile", default="/dev/null")
    parser.add_argument("--n_cvg_bins", default=10, type=int)
    parser.add_argument("--min_CDS", default=50, type=int)
    parser.add_argument("--min_3p_UTR", default=50, type=int)
    parser.add_argument("--min_5p_UTR", default=0, type=int)

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
                                                o.individual_gene,
                                                genes, 
                                                tr_contig,
                                                o.min_CDS,
                                                o.min_3p_UTR,
                                                o.min_5p_UTR)
    
    summary_outrows = []
    summary_simple_outrows = []
    summary_by_exon_outrows = []

    binned_cvg_outrows = []
    CDS_start_stop_outrows = []
    max_t=0
    total_assessed = 0
    for contig, cvg_objs in cvg_objs_by_contig.items():
        print("current contig: %s\tassessed: %d"%(contig, total_assessed))
        for i, cvg_obj in enumerate(cvg_objs):
            total_assessed +=1
            
            t = time.time()
            #if cvg_obj.g.beg != 100706650: continue
            cvg_obj.get_cvg(bamfile)
            if time.time()-t>max_t: max_t = time.time()-t
            #print("t=%f, max_t=%f"%(time.time()-t,max_t))
            #cvg_obj.print_summary()
            summary_outrows.extend(cvg_obj.get_summary_dicts())
            summary_simple_outrows.append(cvg_obj.get_simple_summary_dict())
            summary_by_exon_outrows.extend(cvg_obj.get_by_exon_dicts())

            if o.binned_cvg:
                binned_cvg_outrows.extend(cvg_obj.get_binned_cvg_dicts(o.n_cvg_bins))
            if o.CDS_start_stop:
                CDS_start_stop_outrows.extend(cvg_obj.get_CDS_start_stop_dicts())
            #if i>10:
            #    break
        #break
        
    T_summary = pd.DataFrame(summary_outrows)
    T_summary.to_csv(o.fn_out_summary, index=False, sep="\t")
    
    T_summary_simple = pd.DataFrame(summary_simple_outrows)
    T_summary_simple.to_csv(o.fn_out_summary_simple, index=False, sep="\t")
    
    T_summary_by_exon = pd.DataFrame(summary_by_exon_outrows)
    T_summary_by_exon.to_csv(o.fn_out_summary_by_exon, index=False, sep="\t")

    if o.binned_cvg:
        T_binned_cvg = pd.DataFrame(binned_cvg_outrows)
        T_binned_cvg.to_csv(o.fn_out_binned_cvg, index=False, sep="\t")

    if o.CDS_start_stop:
        T_CDS_start_stop = pd.DataFrame(CDS_start_stop_outrows)
        T_CDS_start_stop.to_csv(o.fn_out_CDS_start_stop, index=False, sep="\t")
        




