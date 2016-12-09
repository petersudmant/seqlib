import argparse
from gtf_to_genes import *
import numpy as np
import scipy.stats as scp_stats
import pandas as pd

import tables

from coveragedata import *
from h5fullgenecvg import h5FullGeneCvg_writer

import logging
import pysam
import pysamstats
import pdb
import math
import time
import sys


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_gtf_index", required=True)
    parser.add_argument("--fn_bam", required=True)
    parser.add_argument("--fn_fasta", required=True)
    parser.add_argument("--meta_type", required=True)
    
    parser.add_argument("--fn_out_summary_simple", required=False)
    parser.add_argument("--fn_out_summary", required=False)
    parser.add_argument("--fn_out_summary_by_exon", required=False)
    parser.add_argument("--fn_out_binned_cvg", required=False)
    parser.add_argument("--fn_out_CDS_start", required=False)
    parser.add_argument("--fn_out_CDS_stop", required=False)
    parser.add_argument("--fn_out_5p", required=False)
    parser.add_argument("--fn_out_3p", required=False)
    parser.add_argument("--fn_out_full_h5", required=False)
    
    parser.add_argument("--individual_gene", required=False, default=None)
    parser.add_argument("--gtf_ID", required=True)
    parser.add_argument("--contig_subset", default=None)
    parser.add_argument("--fn_logfile", default="/dev/null")
    parser.add_argument("--n_cvg_bins", default=10, type=int)
    
    o = parser.parse_args()
    
    logger = logging.getLogger(o.fn_logfile)

    bamfile = pysam.AlignmentFile(o.fn_bam, 'rb')
    
    contig_lengths = {x[0]:x[1] for x in zip(bamfile.references, bamfile.lengths)}
    
    species_id, gtf_path, genes = get_indexed_genes_for_identifier(o.fn_gtf_index,
                                                                   logger, 
                                                                   o.gtf_ID)
    contig_subset = None
    if o.contig_subset: 
        contig_subset = o.contig_subset.split(":")
    else:
        contig_subset = ["chr%d"%i for i in range(1,24)]
    
    cvg_objs_by_contig = get_cvg_objs_by_contig(genes,
                                                "constitutive_single_stop")
    summary_outrows = []
    summary_simple_outrows = []
    summary_by_exon_outrows = []

    binned_cvg_outrows = []
    
    CDS_start_outrows = []
    CDS_stop_outrows = []
    
    end_3p_outrows = []
    end_5p_outrows = []
    
    full_h5 = None
    
    if o.fn_out_full_h5:
        full_h5 = h5FullGeneCvg_writer(o.fn_out_full_h5)

    max_t=0
    total_assessed = 0
    for contig, cvg_objs in cvg_objs_by_contig.items():
        print("current contig: %s\tassessed: %d"%(contig, total_assessed))
        t=time.time()
        """
        ###
        UPDATE WHEN FIXED ON GITHUB 
        ###
        """
        cvg_recarray = pysamstats.load_nondel_coverage(bamfile, 
                                                      chrom=contig, 
                                                      start=0, 
                                                      end=contig_lengths[contig])
        print("time to load contig: %fs"%(time.time()-t))
        sys.stdout.flush()

        for i, cvg_obj in enumerate(cvg_objs):
            total_assessed +=1
            
            t = time.time()
            cvg_obj.get_cvg(cvg_recarray, bamfile)
            if time.time()-t>max_t: 
                max_t = time.time()-t
                print("t=%f, max_t=%f"%(time.time()-t,max_t))
                sys.stdout.flush()
            
            summary_simple_outrows.append(cvg_obj.get_simple_summary_dict())
           

            if o.fn_out_full_h5:
                full_h5.extend(cvg_obj)
            if o.fn_out_summary: 
                summary_outrows.extend(cvg_obj.get_summary_dicts())
            if o.fn_out_summary_by_exon: 
                summary_by_exon_outrows.extend(cvg_obj.get_by_exon_dicts())
            if o.fn_out_binned_cvg:
                binned_cvg_outrows.extend(cvg_obj.get_binned_cvg_dicts(o.n_cvg_bins))
            if o.fn_out_CDS_start:
                CDS_start_outrows.extend(cvg_obj.get_CDS_start_dicts())
            if o.fn_out_CDS_stop:
                CDS_stop_outrows.extend(cvg_obj.get_CDS_stop_dicts())
            if o.fn_out_5p:
                end_5p_outrows.extend(cvg_obj.get_5p_dicts())
            if o.fn_out_3p:
                end_3p_outrows.extend(cvg_obj.get_3p_dicts())
        
    if o.fn_out_summary_simple: 
        T_summary_simple = pd.DataFrame(summary_simple_outrows)
        T_summary_simple.to_csv(o.fn_out_summary_simple, index=False, sep="\t")
    
    if o.fn_out_summary: 
        T_summary = pd.DataFrame(summary_outrows)
        T_summary.to_csv(o.fn_out_summary, index=False, sep="\t")
    
    if o.fn_out_summary_by_exon: 
        T_summary_by_exon = pd.DataFrame(summary_by_exon_outrows)
        T_summary_by_exon.to_csv(o.fn_out_summary_by_exon, index=False, sep="\t")

    if o.fn_out_binned_cvg:
        T_binned_cvg = pd.DataFrame(binned_cvg_outrows)
        T_binned_cvg.to_csv(o.fn_out_binned_cvg, index=False, sep="\t")

    if o.fn_out_CDS_start:
        T_CDS_start = pd.DataFrame(CDS_start_outrows)
        T_CDS_start.to_csv(o.fn_out_CDS_start, index=False, sep="\t")
    
    if o.fn_out_CDS_stop:
        T_CDS_stop = pd.DataFrame(CDS_stop_outrows)
        T_CDS_stop.to_csv(o.fn_out_CDS_stop, index=False, sep="\t")
        
    if o.fn_out_5p:
        T_5p = pd.DataFrame(end_5p_outrows)
        T_5p.to_csv(o.fn_out_5p, index=False, sep="\t")
    
    if o.fn_out_3p:
        T_3p = pd.DataFrame(end_3p_outrows)
        T_3p.to_csv(o.fn_out_3p, index=False, sep="\t")

    if o.fn_out_full_h5:
        full_h5.close()

