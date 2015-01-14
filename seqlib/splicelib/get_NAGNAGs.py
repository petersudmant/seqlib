
"""
GET the NAGNAGs from a junction set
"""

import argparse
import pandas as pd
import pdb
from sys import stderr
import numpy as np


def get_nagnags(all_juncs, strand, FOUT_juncs, FOUT_inf):
    
    if strand == "-": 
        ss_5p_key = "j_right" 
        ss_3p_key = "j_left"
    else:
        ss_5p_key = "j_left" 
        ss_3p_key = "j_right"

    juncs = all_juncs[all_juncs["strand"]==strand]
    juncs = juncs.sort(['contig',ss_5p_key, ss_3p_key])
    prev_ss_5p = None
    prev_ss_3p = None
    for i, row in juncs.iterrows():
        ss_5p = row[ss_5p_key]
        ss_3p = row[ss_3p_key]
        contig = row['contig']
        if prev_ss_5p and ss_5p == prev_ss_5p:
            #print contig, ss_5p, ss_3p, prev_ss_3p, 
            if strand == "+":
                j_left1, j_right1 = ss_5p, prev_ss_3p
                j_left2, j_right2 = ss_5p, ss_3p
            else:
                j_left1, j_right1 = prev_ss_3p, ss_5p
                j_left2, j_right2 = ss_3p, ss_5p 
            size = abs(ss_3p-prev_ss_3p)
            
            FOUT_inf.write("{contig}\t{ss_5p}\t{ss_3p_1}\t{ss_3p_2}\t{strand}\t{size}\n".format(contig=contig,
                                                                                                ss_5p=ss_5p,
                                                                                                ss_3p_1=ss_3p,
                                                                                                ss_3p_2=prev_ss_3p,
                                                                                                strand=strand,
                                                                                                size=size))
            FOUT_juncs.write("{contig}\t{j_left}\t{j_right}\t{strand}\n".format(contig=contig,
                                                                      j_left=j_left1,
                                                                      j_right=j_right1,
                                                                      strand=strand))
            FOUT_juncs.write("{contig}\t{j_left}\t{j_right}\t{strand}\n".format(contig=contig,
                                                                      j_left=j_left2,
                                                                      j_right=j_right2,
                                                                      strand=strand))
            prev_ss_5p = ss_5p
            prev_ss_3p = ss_3p
        else: 
            prev_ss_5p = ss_5p
            prev_ss_3p = ss_3p


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_out_juncs")
    parser.add_argument("--fn_out_info")
    parser.add_argument("--fn_input_juncs")
    parser.add_argument("--contig", default=None)

    o = parser.parse_args()

    juncs = pd.read_csv(o.fn_input_juncs, 
                        header=None, 
                        delimiter="\t",
                        names=["contig",
                               "j_left",
                               "j_right",
                               "strand"])
    FOUT_juncs = open(o.fn_out_juncs,'w') 
    FOUT_inf = open(o.fn_out_info,'w') 
    FOUT_inf.write("contig\tss_5p\tss_3p_1\tss_3p_2\tstrand\tsize\n")

    if o.contig:
        juncs = juncs[juncs["contig"]==o.contig]
        
    get_nagnags(juncs, "-", FOUT_juncs, FOUT_inf)
    get_nagnags(juncs, "+", FOUT_juncs, FOUT_inf)

    

