import argparse
import pdb
import time
import numpy as np
from fastahack import FastaHack
from collections import defaultdict

import splicelib.splicegraph as sg

from junction_writer import JunctionWriter

#import re
#import operator

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_fasta")
    parser.add_argument("--fn_input_gff")
    parser.add_argument("--fn_output_dir", default="output")
    o = parser.parse_args()
    
    fa = FastaHack(o.fn_fasta)
    
    contigs = [contig for contig in fa.names]
    contigs = ["chr2"]
    FOUT = open("%s/exon_sizes.dataframe"%o.fn_output_dir,'w')  
    FOUT_2 = open("%s/exon_sizes_simple.dataframe"%o.fn_output_dir,'w')  
    splice_graphs_by_contig = sg.init_splice_graphs_from_gff3(o.fn_input_gff, contigs=contigs)
    
    FOUT.write("size\n")
    FOUT_2.write("size\n")
    for contig, splice_graph in splice_graphs_by_contig.iteritems():
        for e, v in splice_graph.all_uniq_exons.iteritems():
            FOUT_2.write("%d\n"%(e[1]-e[0]))
        
        for s, es in splice_graph.F_exon_s_e.iteritems():
            p=s
            for e in es:
                if not tuple([s,e]) in splice_graph.all_uniq_exons:
                    print s,e
                    pdb.set_trace()
                FOUT.write("%d\n"%(e-s))
                p=e
        for s, es in splice_graph.R_exon_s_e.iteritems():
            p=s
            for e in es:
                FOUT.write("%d\n"%(e-s))
                p=e
