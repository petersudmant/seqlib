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
    parser.add_argument("--fn_output_gff")
    parser.add_argument("--fn_output_novel_junc_gff")
    parser.add_argument("--fn_output_bed")
    parser.add_argument("--fn_output_juncs")
    parser.add_argument("--track_desc")
    o = parser.parse_args()
    
    fa = FastaHack(o.fn_fasta)
    
    contigs = [contig for contig in fa.names]
    #contigs = ["chr20"]
    print contigs
    
    splice_graphs_by_contig = sg.init_splice_graphs_from_gff3(o.fn_input_gff, contigs=contigs)
        
    F_gff = open(o.fn_output_gff,'w')
    F_novel_gff = open(o.fn_output_novel_junc_gff,'w')
    F_bed = open(o.fn_output_bed,'w')
    j_writer = JunctionWriter(o.fn_output_juncs)

    F_bed.write("""track name=NAGNNAG description="%s NAGNNAGs" visibility=2\n"""%(o.track_desc))
    
    for contig, splice_graph in splice_graphs_by_contig.iteritems():
        seq = fa.get_sequence(contig)
        splice_graph.enumerate_splice_junctions(seq)
        splice_graph.get_NAGNAGs(seq, F_gff, F_novel_gff, F_bed, j_writer)
    

