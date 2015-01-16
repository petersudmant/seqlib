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
    parser.add_argument("--fn_output_gff_dir", default="gff")
    parser.add_argument("--fn_output_bed_dir", default="bed")
    parser.add_argument("--fn_output_juncs_dir", default="junctions")
    parser.add_argument("--track_desc")
    parser.add_argument("--micro_exon_type", choices=['skipped','5p','3p'], required=True)
    parser.add_argument("--max_alt_exon_len", default=60, type=int)
    o = parser.parse_args()
    
    fa = FastaHack(o.fn_fasta)
    
    contigs = [contig for contig in fa.names]
    contigs = ["chr2"]
    print contigs
    
    splice_graphs_by_contig = sg.init_splice_graphs_from_gff3(o.fn_input_gff, contigs=contigs)
    
    mx_type = o.micro_exon_type 
    F_gff = open("{gff_dir}/{mx_type}.gff".format(gff_dir=o.fn_output_gff_dir,
                                                  mx_type=mx_type),'w')
    F_novel_gff = open("{gff_dir}/{mx_type}.novel.gff".format(gff_dir=o.fn_output_gff_dir,
                                                              mx_type=mx_type),'w')
    F_bed = open("{bed_dir}/{mx_type}.bed".format(bed_dir=o.fn_output_bed_dir,
                                                        mx_type=mx_type),'w')
    track_desc="""track name="{mx_type}" description="{track_desc} {mx_type}" visibility=2 itemRgb=On\n""".format(track_desc=o.track_desc,
                                                                                                  mx_type=mx_type)
    F_bed.write(track_desc)

    j_writer = JunctionWriter("{junction_dir}/{mx_type}".format(junction_dir=o.fn_output_juncs_dir,
                                                                 mx_type=mx_type))

    for contig, splice_graph in splice_graphs_by_contig.iteritems():
        seq = fa.get_sequence(contig)
        #splice_graph.enumerate_splice_junctions(seq)
        if mx_type=="5p" or mx_type=="3p":
            get_3p=mx_type=="3p"
            get_5p=mx_type=="5p"
            splice_graph.get_5p_3p_alt_exon(seq, 
                                            F_gff, 
                                            F_novel_gff, 
                                            F_bed, 
                                            j_writer, 
                                            get_5p=get_5p,
                                            get_3p=get_3p,
                                            n=o.max_alt_exon_len)
        else:
            print "NOT IMPLEMENTED"

