import argparse
import pdb
import time
import numpy as np
from fastahack import FastaHack
from collections import defaultdict
import BCBio.GFF as GFF
import collections

from bx.intervals.cluster import ClusterTree
import time
from junction_writer import JunctionWriter
import pysam
import pandas as pd
from sys import stderr

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_out")
    parser.add_argument("--fn_filtered_juncs")
    parser.add_argument("--min_entropy", default=1.5, type=float)
    parser.add_argument("--fn_inputs", nargs="*")
    o = parser.parse_args()
   
    tables_by_tissue = {}
    pileup_tables_by_tissue = {}
    entropy_cols = []
    overhang_cols = []
    pileup_cols = []
    for f in o.fn_inputs:
        tissue = f.split("/")[-1].split(".")[0]
        t = pd.read_csv(f, 
                        header=None, 
                        delimiter="\t",
                        names=["contig",
                               "j_left",
                               "j_right",
                               "strand",
                               "%s_entropy"%tissue,
                               "%s_min_overhang"%tissue])
        t_pileup = pd.read_csv("%s.pileup"%f, 
                                header=None, 
                                skiprows=1, 
                                delimiter="\t",
                                names=["position", "read_count_%s"%tissue])
        pileup_tables_by_tissue[tissue] = t_pileup
        pileup_cols.append("read_count_%s"%tissue)
        entropy_cols.append("%s_entropy"%tissue)
        overhang_cols.append("%s_min_overhang"%tissue)
        tables_by_tissue[tissue] = t
        print >> stderr, tissue   
    
    tables = tables_by_tissue.values() 
    T = tables[0]
    for i in xrange(1,len(tables)):
        t = tables[i]
        T = pd.merge(T, t, on=["contig", "j_left","j_right","strand"])
    
    pileup_tables = pileup_tables_by_tissue.values() 
    pileup_T = pileup_tables[0]
    for i in xrange(1,len(pileup_tables)):
        t = pileup_tables[i]
        pileup_T = pd.merge(pileup_T, t, on=["position"],how='outer')
    pielup_T = pileup_T.fillna(0)
    pileup_T.to_csv("%s.pileups"%o.fn_out, sep="\t", index=False)

    cols = T.columns
    e_cols = sorted([np.where(cols==e)[0][0] for e in entropy_cols])
    overhang_cols = sorted([np.where(cols==e)[0][0] for e in overhang_cols])
    
    all_ents = T.values[:,e_cols]
    all_overhangs = T.values[:,overhang_cols]

    Tout=T[~(np.sum(all_ents==0,1)==len(e_cols))]
    #concat_T = pd.concat(tables_by_tissue.values())
    #concat_T.to_csv(o.fn_out, sep="\t", index=False)
    Tout.to_csv(o.fn_out, sep="\t", index=False)
    
    T_switch = T[np.sum(all_ents<o.min_entropy,1)==(len(e_cols)-1)]
    T_switch.to_csv("%s.one_gt_min_ent"%o.fn_out, sep="\t", index=False)

    T_filtered = T[(np.sum(all_ents>=o.min_entropy,1)>=1)]
    T_filtered.to_csv(o.fn_filtered_juncs, 
                      sep="\t", 
                      index=False, 
                      columns = ["contig", "j_left", "j_right", "strand"], 
                      header=False) 


      
    #OUTPUT ALL JUNCTIONS WITH at least one min_overhang >6 AND one entropy > 2
