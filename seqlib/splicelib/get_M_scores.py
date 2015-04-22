import argparse
import pdb
import os
import numpy as np
import pandas as pd
from collections import defaultdict


def get_M_scores(sample_order, t_by_comparison):

    print(sample_order)
    pdb.set_trace()

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--BF_cutoff", type=int, default=10)
    parser.add_argument("--fn_MISO_comparisons", nargs="+")
    parser.add_argument("--sample_order", nargs="+")
    o = parser.parse_args()
     
    for fn_compare in o.fn_MISO_comparisons:
        comp = fn_compare.split("/")[-2]
        s1, s2 = comp.split("_vs_")
        t=pd.read_csv(fn_compare, header=0, sep="\t")
        t = t[['event_name','bayes_factor','sample1_posterior_mean','sample2_posterior_mean']]
        t=t.set_index(['event_name'])
        t_by_comparison[tuple([s1,s2])] = t
    
    Ms = get_M_scores(o.sample_order, t_by_comparison)
