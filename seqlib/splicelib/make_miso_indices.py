import argparse
import pdb
import time
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
    o = parser.parse_args()
    
    fa = FastaHack(o.fn_fasta)
    contigs = [contig for contig in fa.names]

    contigs = ["chr19"]
    if o.force_contig:
        contigs = [c for c in contigs if c==o.force_contig]
    print contigs
    
    contig_sizes = {}
    for contig in contigs:
        contig_sizes[contig] = fa.get_sequence_length(contig)

    splice_graphs_by_contig = sg.init_splice_graphs_from_gff3(o.fn_input_gff, 
                                                              contigs=contigs)
    m_util = mu.MisoUtils(sg_by_contig = splice_graphs_by_contig, contig_sizes = contig_sizes)

    m_util.define_SE_events("{outdir}/SE.gff".format(outdir=o.fn_out_dir))
    m_util.define_A3SS_events("{outdir}/A3SS.gff".format(outdir=o.fn_out_dir))
    m_util.define_A5SS_events("{outdir}/A5SS.gff".format(outdir=o.fn_out_dir))
    m_util.define_MXE_events("{outdir}/MXE.gff".format(outdir=o.fn_out_dir))
    m_util.define_RI_events("{outdir}/RI.gff".format(outdir=o.fn_out_dir))
    m_util.define_ALE_events("{outdir}/ALE.gff".format(outdir=o.fn_out_dir))
    m_util.define_AFE_events("{outdir}/AFE.gff".format(outdir=o.fn_out_dir))
    
