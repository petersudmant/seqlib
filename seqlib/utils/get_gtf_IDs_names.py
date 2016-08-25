import argparse
from gtf_to_genes import *
import numpy as np
import scipy.stats as scp_stats
import pandas as pd

import tables

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
    parser.add_argument("--gtf_ID", required=True)
    parser.add_argument("--fn_logfile", default="/dev/null")
    parser.add_argument("--fn_out")
    
    o = parser.parse_args()
    
    logger = logging.getLogger(o.fn_logfile)
    
    species_id, gtf_path, genes = get_indexed_genes_for_identifier(o.fn_gtf_index,
                                                                   logger, 
                                                                   o.gtf_ID)
    
    outrows = []
    for g in genes['protein_coding']:
        name = g.names[0]
        gene_id = g.gene_id
        for t in g.transcripts:
            outrows.append({"gene_id":gene_id,
                            "transcript_id":t.cdna_id,
                            "name":name})
    
    
    t = pd.DataFrame(outrows)
    t.to_csv(o.fn_out, sep="\t", index=False)
