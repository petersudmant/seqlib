import argparse
import pdb
import os
import numpy as np
from fastahack import FastaHack
from collections import defaultdict

import splicelib.splicegraph as sg
import splicelib.misoutils as mu

from junction_writer import JunctionWriter
#import re
#import operator

if __name__=="__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--fn_fasta")
    parser.add_argument("--fn_input_gff")
    parser.add_argument("--fn_out_dir", default="miso_index")
    parser.add_argument("--force_contig", default=None)
    parser.add_argument("--force_feature", default=None)
    o = parser.parse_args()
    
    fa = FastaHack(o.fn_fasta)
    contigs = [contig for contig in fa.names]
    
    if o.force_contig:
        contigs = [c for c in contigs if c==o.force_contig]
    print contigs
    
    contig_sizes = {}
    for contig in contigs:
        contig_sizes[contig] = fa.get_sequence_length(contig)
    
    if not os.path.exists("{outdir}/gff".format(outdir=o.fn_out_dir)):
        os.mkdir("{outdir}/gff".format(outdir=o.fn_out_dir))
    if not os.path.exists("{outdir}/bed".format(outdir=o.fn_out_dir)):
        os.mkdir("{outdir}/bed".format(outdir=o.fn_out_dir))
    if not os.path.exists("{outdir}/info".format(outdir=o.fn_out_dir)):
        os.mkdir("{outdir}/info".format(outdir=o.fn_out_dir))
    
    if o.force_feature:
        splice_graphs_by_contig = sg.init_splice_graphs_from_gff3(o.fn_input_gff, 
                                                                  contigs=contigs,
                                                                  features=[o.force_feature])
    else:
        splice_graphs_by_contig = sg.init_splice_graphs_from_gff3(o.fn_input_gff, 
                                                                  contigs=contigs)


    m_util = mu.MisoUtils(sg_by_contig = splice_graphs_by_contig,
                          contig_sizes = contig_sizes,
                          fn_info = "{outdir}/info/info.df".format(outdir=o.fn_out_dir))

    m_util.define_SE_events("{outdir}/gff/SE.gff".format(outdir=o.fn_out_dir),
                            "{outdir}/bed/SE.bed".format(outdir=o.fn_out_dir))
    
    m_util.define_A3SS_events("{outdir}/gff/A3SS.gff".format(outdir=o.fn_out_dir),
                              "{outdir}/bed/A3SS.bed".format(outdir=o.fn_out_dir))
    
    m_util.define_A5SS_events("{outdir}/gff/A5SS.gff".format(outdir=o.fn_out_dir),
                              "{outdir}/bed/A5SS.bed".format(outdir=o.fn_out_dir))
    
    m_util.define_MXE_events("{outdir}/gff/MXE.gff".format(outdir=o.fn_out_dir),
                              "{outdir}/bed/MXE.bed".format(outdir=o.fn_out_dir))
    
    m_util.define_AFE_events("{outdir}/gff/AFE.gff".format(outdir=o.fn_out_dir),
                              "{outdir}/bed/AFE.bed".format(outdir=o.fn_out_dir))
    
    m_util.define_ALE_events("{outdir}/gff/ALE.gff".format(outdir=o.fn_out_dir),
                             "{outdir}/bed/ALE.bed".format(outdir=o.fn_out_dir))
    
    m_util.define_RI_events("{outdir}/gff/RI.gff".format(outdir=o.fn_out_dir),
                            "{outdir}/bed/RI.bed".format(outdir=o.fn_out_dir))
    
    m_util.output_info()
