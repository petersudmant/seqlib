import tables
import logging 
import numpy as np
import pandas as pd
import argparse
import sys
import pdb
import pysam
import scipy.stats as stats

import itertools

from expression.coveragedata import *
from gtf_to_genes import *
import timeit



if __name__=="__main__":

    parser = argparse.ArgumentParser()
    
    #output feature centered counts
    parser.add_argument("--fn_out", required=True)
    parser.add_argument("--fn_gtf_index", required=True)
    parser.add_argument("--gtf_ID", required=True)
    parser.add_argument("--fn_logfile", default='/dev/stderr')
    parser.add_argument("--contig", default=None)
    
    args = parser.parse_args()
    
    sys.stderr.write("loading gene annotations...")
    logger = logging.getLogger(args.fn_logfile)
    s_id, path, genes = get_indexed_genes_for_identifier(args.fn_gtf_index,
                                                         logger, 
                                                         args.gtf_ID)
    sys.stderr.write("done\n")
    cvg_objs_by_contig = get_cvg_objs_by_contig(genes,
                                                "transcript", 
                                                contig_subset=args.contig)
    
    starts = []
    ends = []
    contigs = []
    strands = []
    types = []
    
    genes = []
    transcripts = []
    attributes = []

    outrows = []

    attribute_str = """transcript_id "{TID}"; gene_id "{gene_id}"; gene_name "{gene_name}";"""
    for contig, cvg_obs in cvg_objs_by_contig.items():
        for cvg_ob in cvg_obs:
            for ex in cvg_ob.exons: 
                starts.append(ex[0])
                ends.append(ex[1])
                contigs.append(contig)
                strands.append(cvg_ob.strand and "+" or "-")
                genes.append(cvg_ob.gene_id)
                transcripts.append(cvg_ob.TID)
                types.append("exon")
                attributes.append(attribute_str.format(TID=cvg_ob.TID,
                                                      gene_id=cvg_ob.gene_id,
                                                      gene_name=cvg_ob.g.names[0]))

            for intr in get_introns(cvg_ob.exons): 
                starts.append(intr[0])
                ends.append(intr[1])
                contigs.append(contig)
                strands.append(cvg_ob.strand and "+" or "-")
                genes.append(cvg_ob.gene_id)
                transcripts.append(cvg_ob.TID)
                types.append("intron")
                attributes.append(attribute_str.format(TID=cvg_ob.TID,
                                                      gene_id=cvg_ob.gene_id,
                                                      gene_name=cvg_ob.g.names[0]))
    
    T = pd.DataFrame({"contig":contig, 
                      "source": types,
                      "feature": types,
                      "start":starts,
                      "end": ends,
                      "strand":strands,
                      "attribute":attributes})
    T['score'] = "."
    T['frame'] = "."
    columns =  ["contig", 
                "source", 
                "feature", 
                "start",
                "end",
                "score", 
                "strand", 
                "frame", 
                "attribute"]

    T.to_csv(args.fn_out, sep="\t", columns=columns, header=False, index=False)


