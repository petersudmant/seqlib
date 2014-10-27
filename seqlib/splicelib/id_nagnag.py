import argparse
import pdb
import time
import numpy as np
from fastahack import FastaHack
from collections import defaultdict

import splicelib.splicegraph as sg

#import re
#import operator

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_fasta")
    parser.add_argument("--fn_gff")
    o = parser.parse_args()
    
    fa = FastaHack(o.fn_fasta)
    
    contigs = [contig for contig in fa.names]
    contigs = ["19"]
    #contigs = ["GL000213.1"]
    splice_graphs_by_contig = sg.init_splice_graphs_from_gff(o.fn_gff, contigs)

    for contig, splice_graph in splice_graphs_by_contig.iteritems():
        seq = fa.get_sequence(contig)
        splice_graph.enumerate_splice_junctions(seq)
    

