import argparse
from gtf_to_genes import *
import numpy as np
import scipy.stats as scp_stats
import pandas as pd
import pdb

import logging
import pysam
import pysamstats
import pdb
import math
import time
import sys

import scipy.stats as scp_stats

def get_idxstats(fn_bam):
    n_mapped_by_contig = {}
    n_mapped_total = 0
    for inf in pysam.idxstats(fn_bam):
        contig,len,mapped,unmapped = inf.split()
        mapped=int(mapped)
        n_mapped_by_contig[contig] = mapped
        n_mapped_total +=mapped
    
    return n_mapped_total, n_mapped_by_contig

def running_mean(x, N):
    cumsum = np.cumsum(np.insert(x, 0, 0))
    return (cumsum[N:] - cumsum[:-N]) / N 

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_bam", required=True)
    parser.add_argument("--locus", default=None)
    parser.add_argument("--fn_out", default=None)
    parser.add_argument("--tile", default=None, type=int)

    o = parser.parse_args()
    
    bamfile = pysam.AlignmentFile(o.fn_bam, 'rb')
    
    n_mapped_total, n_mapped_by_contig = get_idxstats(o.fn_bam)

    contig, s_e = o.locus.split(":")
    start, end = s_e.split("-")
    start, end = int(start), int(end)
    
    ## #NOTE - still has maxdepth issue
    max_depth=32000
    #max_depth=64000
    #max_depth=256000
    while(1):
        print("step", max_depth)
        cvg_recarray = pysamstats.load_nondel_coverage(bamfile, 
                                                      chrom=contig, 
                                                      start=start, 
                                                      end=end, 
                                                      pad=True,
                                                      max_depth=max_depth)
        print(np.amax(cvg_recarray.reads_all))
        
        if np.amax(cvg_recarray.reads_all)>=max_depth-max_depth/3.0:
            max_depth = max_depth*2
        else:
            break
    
    mu_cvg = np.mean(cvg_recarray.reads_all)
    
    outrows = []
    if o.tile:
        #smoothed = pd.rolling_mean(cvg_recarray.reads_all,o.tile)
        smoothed = cvg_recarray.reads_all
        for i in range(start, end, o.tile):
            offset=i*o.tile
            outrows.append({"n_mapped_total":n_mapped_total,
                            "n_mapped_to_contig":n_mapped_by_contig[contig],
                            "contig":contig,
                            "start":start+offset,
                            "end":start+offset+o.tile,
                            "cvg":smoothed[i],
                            "mu_cvg": mu_cvg})
    else:
        outrows.append({"n_mapped_total":n_mapped_total,
                        "n_mapped_to_contig":n_mapped_by_contig[contig],
                        "contig":contig,
                        "start":start,
                        "end":end,
                        "cvg": mu_cvg})
    
    t = pd.DataFrame(outrows) 
    t.to_csv(o.fn_out, sep="\t", index=False)

