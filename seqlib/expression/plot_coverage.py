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

    
if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_gtf_index", required=True)
    parser.add_argument("--fn_bams", required=True, nargs="+")
    parser.add_argument("--samples", required=True, nargs="+")
    parser.add_argument("--fn_out", required=True)
    parser.add_argument("--gene", required=True, default=None)
    parser.add_argument("--gtf_ID", required=True)
    parser.add_argument("--fn_logfile", default="/dev/null")

    o = parser.parse_args()

    logger = logging.getLogger(o.fn_logfile)
    
    bamfiles = {}
    for i, sample in enumerate(o.samples):
        bamfiles[sample] = pysam.AlignmentFile(o.fn_bams[i], 'rb')
        
    species_id, gtf_path, genes = get_indexed_genes_for_identifier(o.fn_gtf_index,
                                                                   logger, 
                                                                   o.gtf_ID)
    
    g_obj = None
    for g in genes['protein_coding']:
        if o.gene in g.names:
            g_obj = g
            break
    
    assert g_obj!=None, "gene name {name} not found".format(o.gene)
    
    bp_cov_tables = []
    for sample in o.samples:
        cvg_obj = CoverageData(g_obj)
        bp_cvg_rows = cvg_obj.get_cvg(bamfiles[sample])
        bp_cvg_rows = cvg_obj.get_bp_cvg_dicts()
        T_bp_cov = pd.DataFrame(bp_cvg_rows)
        T['sample'] = sample
        bp_cov_tables.append(T)
    
    T = pd.concat(bp_cov_tables)
    T.to_csv(o.fn_out, index=False, sep="\t")
