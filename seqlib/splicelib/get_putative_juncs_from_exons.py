import argparse
import pdb
import time
import numpy as np
from fastahack import FastaHack
from collections import defaultdict
import pandas as pd

import splicelib.splicegraph as sg
from junction_writer import JunctionWriter

#import re
#import operator

def get_juncs(contig, s, e, strand, LR_intron_interval_tree):
    
    juncs = []
    for interval in LR_intron_interval_tree.find(s,e):
        if interval.start < s and interval.end > e:
            juncs.append({"contig":contig, "start": interval.start-1, "end":s, "strand":strand})
            juncs.append({"contig":contig, "start": e-1, "end":interval.end, "strand":strand})
            juncs.append({"contig":contig, "start": interval.start-1, "end":interval.end, "strand":strand})
    
    return juncs

if __name__=="__main__":
    """
    outputs tophat formatted exon juncs
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_fasta")
    parser.add_argument("--fn_input_gff")
    parser.add_argument("--fn_input_exons")
    parser.add_argument("--fn_output_juncs")
    parser.add_argument("--force_contig", default=None)
    o = parser.parse_args()
    
    fa = FastaHack(o.fn_fasta)
    
    contigs = [contig for contig in fa.names]
    if o.force_contig:
        contigs = [c for c in contigs if c==o.force_contig]
    print contigs
    
    exons = pd.read_csv(o.fn_input_exons, 
                        header=None, 
                        index_col=None,
                        sep="\t",
                        names=["contig", "start", "end", "strand"])

    splice_graphs_by_contig = sg.init_splice_graphs_from_gff3(o.fn_input_gff, contigs=contigs)
    
    juncs = []
    for contig, splice_graph in splice_graphs_by_contig.iteritems():
        seq = fa.get_sequence(contig)
        #splice_graph.enumerate_splice_junctions(seq)
        curr_exons = exons[exons['contig']==contig]
        
        for i, e_row in curr_exons.iterrows():
            contig, s, e, strand = e_row['contig'], e_row['start'], e_row['end'], e_row['strand']
            if strand == "+":
                juncs += get_juncs(contig, s, e, strand, splice_graph.F_LR_intron_interval_tree)
            else:
                juncs += get_juncs(contig, s, e, strand, splice_graph.R_LR_intron_interval_tree)
    
    T = pd.DataFrame(juncs)
    T.to_csv(o.fn_output_juncs,
             sep="\t", 
             columns = ["contig", "start", "end", "strand"],
             index=None, 
             header=False)
