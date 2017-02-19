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

if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_out")
    parser.add_argument("--fn_gtf_index")
    parser.add_argument("--gtf_ID")
    parser.add_argument("--fn_logfile", default="/dev/null")
    parser.add_argument("--fn_fasta")
    args = parser.parse_args()
    
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
            UTR_3p_seqs = [fa.fetch(contig,UTR_e[0],UTR_e[1]-1) for UTR_e in cvg_obj.UTR_3p_exons]
            UTR_5p_seqs = [fa.fetch(contig,UTR_e[0],UTR_e[1]-1) for UTR_e in cvg_obj.UTR_5p_exons]
            CDS_seqs = [fa.fetch(contig,CDS_e[0],CDS_e[1]-1) for CDS_e in cvg_obj.coding_exons]
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
            
            outrows.append({"gene_id":g.gene_id,
                            "gene_name":g.names[0],
                            "transcript_id":cvg_obj.TID,
                            "contig":contig,
                            "strand":str_strand,
                            "seq":seq,
                            "CDS_start":len(UTR_5p_seq),
                            "CDS_end":len(UTR_5p_seq)+len(CDS_seq),
                            "n_exons":n_exons})
    
    T = pd.DataFrame(outrows)
    T.to_csv(args.fn_out, sep="\t", index=False,compression="gzip")

