import argparse
import logging
import pandas as pd
from gtf_to_genes import *
from expression.coveragedata import *
import pdb
import string
import re
import pysam
import sys

trans = string.maketrans('ATCGatcg', 'TAGCtacg')
def revcomp(s):
    return s.translate(trans)[::-1]

def get_transcript_seq_inf(fa, contig, cvg_obj):

    UTR_3p_seqs = [fa.fetch(contig,UTR_e[0],UTR_e[1]) for UTR_e in cvg_obj.UTR_3p_exons]
    UTR_5p_seqs = [fa.fetch(contig,UTR_e[0],UTR_e[1]) for UTR_e in cvg_obj.UTR_5p_exons]
    CDS_seqs = [fa.fetch(contig,CDS_e[0],CDS_e[1]) for CDS_e in cvg_obj.coding_exons]
    n_exons = len(cvg_obj.exons) 
    g = cvg_obj.g
    if g.strand: 
        UTR_3p_seq = "".join(UTR_3p_seqs)
        UTR_5p_seq = "".join(UTR_5p_seqs)
        CDS_seq = "".join(CDS_seqs)
        str_strand = "+"
    else:
        UTR_3p_seq = "".join([revcomp(s) for s in UTR_3p_seqs[::-1]])
        UTR_5p_seq = "".join([revcomp(s) for s in UTR_5p_seqs[::-1]])
        CDS_seq = "".join([revcomp(s) for s in CDS_seqs[::-1]])
        str_strand = "-"
    
    seq = UTR_5p_seq + CDS_seq + UTR_3p_seq
        
    return {"gene_id":g.gene_id,
            "gene_name":g.names[0],
            "transcript_id":cvg_obj.TID,
            "contig":contig,
            "strand":str_strand,
            "seq":seq,
            "CDS_seq":CDS_seq,
            "UTR_3p_seq":UTR_3p_seq,
            "UTR_5p_seq":UTR_5p_seq,
            "CDS_start":len(UTR_5p_seq),
            "CDS_end":len(UTR_5p_seq)+len(CDS_seq),
            "n_exons":n_exons}


def get_sequence(args):

    fa = pysam.Fastafile(args.fn_fasta)
    
    sys.stderr.write("loading gene annotations...")
    logger = logging.getLogger(args.fn_logfile)
    s_id, path, genes = get_indexed_genes_for_identifier(args.fn_gtf_index,
                                                         logger, 
                                                         args.gtf_ID)

    sys.stderr.write("done\n")
    cvg_objs_by_contig = get_cvg_objs_by_contig(genes,
                                                "transcript")
    outrows = []
    for contig, cvg_objs in cvg_objs_by_contig.items():
        if not contig in fa.references:
            continue
        sys.stderr.write("{contig}...".format(contig=contig))
        for cvg_obj in cvg_objs:
            seq_inf = get_transcript_seq_inf(fa, contig, cvg_obj)
            outrows.append(seq_inf)
                
    T = pd.DataFrame(outrows)
    T.to_csv(args.fn_out, sep="\t", index=False,compression="gzip")

def get_features(args):
    """
    output all features in genomic coords (g_start, g_end)
    output all features in transcript coords (t_start, t_end)
    
    """

    stop_codons = "(?=TGA|TAA|TGA)"
    
    if args.feature_type == "CDS_CODON":
        codons = "(?=%s)"%("|".join(args.feature_args))

    fa = pysam.Fastafile(args.fn_fasta)
    sys.stderr.write("loading gene annotations...")
    logger = logging.getLogger(args.fn_logfile)
    s_id, path, genes = get_indexed_genes_for_identifier(args.fn_gtf_index,
                                                         logger, 
                                                         args.gtf_ID)
    sys.stderr.write("done\n")
    cvg_objs_by_contig = get_cvg_objs_by_contig(genes,
                                                "transcript")
    outrows = []
    for contig, cvg_objs in cvg_objs_by_contig.items():
        if not contig in fa.references:
            continue
        sys.stderr.write("{contig}...".format(contig=contig))
        for cvg_obj in cvg_objs:
            
            seq_inf = get_transcript_seq_inf(fa, contig, cvg_obj)
            CDS_end = seq_inf['CDS_end']
            CDS_start = seq_inf['CDS_start']
            
            seq = seq_inf['seq'] 
            CDS_seq = seq_inf['CDS_seq']
            UTR_3p_seq = seq_inf['UTR_3p_seq']

            genome_positions = cvg_obj.get_transcript_to_genome_coords()
            
            if args.feature_type == "UTR_STOP":
                codon_ps = [m.start()+CDS_end for m in re.finditer(stop_codons,UTR_3p_seq)]
            elif args.feature_type == "CDS_CODON":
                codon_ps = [m.start()+CDS_start for m in re.finditer(codons,CDS_seq) if m.start() %3==0]
            for i, t_start in enumerate(codon_ps):
                ###NOTE: get the +2 position, then sort, then add 1. Why? because if
                ### on negative strand, then +3 will get the position BEFORE the first base
                ### need to get the 0 based coords, sort them, then add 1 to the last one
                g_start_end = sorted([genome_positions[t_start],genome_positions[t_start+2]])
                g_start = g_start_end[0]
                g_end = g_start_end[1]+1 
                outrows.append({"gene_id": seq_inf['gene_id'],
                                "gene_name": seq_inf['gene_name'],
                                "transcript_id": seq_inf['transcript_id'],
                                "contig": contig,
                                "t_start": t_start,
                                "t_end": t_start+3,
                                "g_start": g_start,
                                "g_end": g_end,
                                "feature_idx":i,
                                "feature_seq":seq[t_start:t_start+3],
                                "strand":cvg_obj.strand})
    T = pd.DataFrame(outrows)
    T.to_csv(args.fn_out, 
             sep="\t", 
             index=False,
             compression="gzip")


if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers() 
    
    #get sequence
    parser_getseq = subparsers.add_parser("getsequence")
    parser_getseq.add_argument("--fn_out", required=True)
    parser_getseq.add_argument("--fn_gtf_index", required=True)
    parser_getseq.add_argument("--gtf_ID", required=True)
    parser_getseq.add_argument("--fn_fasta", required=True)
    parser_getseq.add_argument("--fn_logfile", default="/dev/null")
    parser_getseq.set_defaults(func=get_sequence)
    
    parser_getfeatures = subparsers.add_parser("getfeatures")
    parser_getfeatures.add_argument("--fn_out", required=True)
    parser_getfeatures.add_argument("--fn_gtf_index", required=True)
    parser_getfeatures.add_argument("--gtf_ID", required=True)
    parser_getfeatures.add_argument("--fn_fasta", required=True)
    parser_getfeatures.add_argument("--feature_type", required=True, choices=["UTR_STOP", "CDS_CODON"])
    parser_getfeatures.add_argument("--feature_args", required=False, nargs="+")
    parser_getfeatures.add_argument("--fn_logfile", default="/dev/null")
    parser_getfeatures.set_defaults(func=get_features)

    args = parser.parse_args()
    args.func(args)
    
