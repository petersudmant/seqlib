import numpy as np
import pandas as pd
import argparse
import sys
import pdb
import pysam
    
def get_reads(args):
    
    bamfile = pysam.AlignmentFile(args.fn_bam, 'rb')
    contig = args.contig
    start = args.start
    end = args.end
    
    outrows = []
    for r in bamfile.fetch(contig, start, end):
        strand = not(r.is_reverse)
        pos = r.pos
        qual = r.mapping_quality
        #only take unique reads - STAR -> 255 for uniq
        if qual!=255: continue

        #NOW set pos to 5' end
        if strand:
            pos = pos
        else:
            pos = pos+r.query_alignment_length

        #print(r)
        #print(pos, strand, r.query_alignment_length, r.cigartuples)
        #pdb.set_trace()
        outrows.append({"pos":pos,
                        "strand":strand,
                        "len":r.query_alignment_length,
                        "sample":args.sample})

    t = pd.DataFrame(outrows)
    t.to_csv(args.fn_out, sep="\t", index=False)
    

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers() 
    
    #spit out reads
    parser_create = subparsers.add_parser("getreads")
    parser_create.add_argument("--fn_bam", required=True)
    parser_create.add_argument("--fn_out", required=True)
    parser_create.add_argument("--contig", required=True)
    parser_create.add_argument("--start", required=True, type=int)
    parser_create.add_argument("--end", required=True, type=int)
    parser_create.add_argument("--sample", required=True)
    parser_create.set_defaults(func=get_reads)
    
    args = parser.parse_args()
    args.func(args)
