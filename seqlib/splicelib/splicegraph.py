import numpy as np
import pdb
import argparse
import BCBio.GFF as GFF
import re
import string 

from transcript import Transcript
from collections import defaultdict
from sys import stderr
from fastahack import FastaHack

trans = string.maketrans('ATCGatcg', 'TAGCtacg')

def revcomp(s):
    return s.translate(trans)[::-1]

class SpliceGraph(object):

    def __init__(self, contig, **kwargs):
            
        """
          5'      3'
        + GT------AG
        
          3'      5'
        - CT------AC

        a collection of hashes of lists representing all splice site 
        and all exon boundary relationships

        """
        
        self.contig = contig
        
        self.F_3p_5p_ss = kwargs["F_3p_5p_ss"]
        self.F_5p_3p_ss = kwargs["F_5p_3p_ss"]
        self.R_3p_5p_ss = kwargs["R_3p_5p_ss"]
        self.R_5p_3p_ss = kwargs["R_5p_3p_ss"]
        
        self.F_exon_s_e = kwargs["F_exon_s_e"]
        self.F_exon_e_s = kwargs["F_exon_e_s"]
        self.R_exon_s_e = kwargs["R_exon_s_e"]
        self.R_exon_e_s = kwargs["R_exon_e_s"]
        
        self.F_5p_to_gene_info = kwargs["F_5p_to_gene_info"]
        self.F_3p_to_gene_info = kwargs["F_3p_to_gene_info"]
        self.R_5p_to_gene_info = kwargs["R_5p_to_gene_info"]
        self.R_3p_to_gene_info = kwargs["R_3p_to_gene_info"]
        
    def get_common_shortest_exon(self, ss, ss_type, strand):
        """
        ss_type = 5' | 3'
        strand = FWD | REV 

        all exons should be returned s<e

        NOTE: given 
                5'-1
            -----  5'-2
            --------
            this will return 5'-1
            BUT, given
            
                5'-1
            -----  5'-2
            --------
              -- 5'-3
            this will still return 5'-1, I THINK this is correct...
        """
        assert ss_type == "5'" or ss_type == "3'"
        assert strand == "FWD" or strand == "REV"

        if strand == "FWD":
            if ss_type == "5'":
                return tuple([max(self.F_exon_e_s[ss]), ss])
            else:
                return tuple([ss, min(self.F_exon_s_e[ss])])
        else:
            if ss_type == "5'":
                return tuple([ss, min(self.R_exon_s_e[ss])])
            else:
                return tuple([max(self.R_exon_e_s[ss]), ss])
    

    def get_NAGNAGs(self, seq, F_gff, F_bed, n=33):
        """
        define NAGNAG
        SD SA1,SA2 where SA1 and SA2 are within n bp of each other
        
        assuming common shortest here in all cases
        """
        
        #FWD
        """
        get all NAGNAGs on FWD strand
        """

        for strand_d, _5p_3p_ss in { "FWD" : self.F_5p_3p_ss, "REV" : self.R_5p_3p_ss }.iteritems():
            for ss_5p, ss_3p_list in _5p_3p_ss.iteritems():
                annotated_ss_3p = np.sort(np.array(ss_3p_list))

                max_ss_3p = annotated_ss_3p[-1]
                min_ss_3p = annotated_ss_3p[0]
                
                if strand_d == "FWD":
                    annot_3p_ss = max_ss_3p
                else:
                    annot_3p_ss = min_ss_3p

                ds_exon = self.get_common_shortest_exon(annot_3p_ss, "3'", strand_d)
                us_exon = self.get_common_shortest_exon(ss_5p, "5'", strand_d)
                
                if strand_d == "FWD":
                    us_gene_inf = self.F_5p_to_gene_info[ss_5p]
                    ds_gene_inf = self.F_3p_to_gene_info[annot_3p_ss]
                else:
                    us_gene_inf = self.R_5p_to_gene_info[ss_5p]
                    ds_gene_inf = self.R_3p_to_gene_info[annot_3p_ss]

                gene_inf = []
                for inf in us_gene_inf+ds_gene_inf:
                    if not inf in gene_inf: gene_inf.append(inf)
            
                """
                k is the number of bp up/down of the annotated 
                3p ss to search for now ones
                """
                
                k = min(n, abs(annot_3p_ss-ss_5p))
                ss_seq = strand_d == "FWD" and "AG" or "CT"
                if strand_d == "FWD":
                    seq_s, seq_e = annot_3p_ss-k, annot_3p_ss-2
                else:
                    seq_s, seq_e = annot_3p_ss+2, annot_3p_ss+k
                """
                offset by a 2 if in fwd dir as you want exon 
                EXON start, not the, position of the ss
                for rev, the position is fine already
                """
                delta = strand_d == "FWD" and 2 or 0 
                
                alt_3p_ss = np.array([m.start()+seq_s+delta for m in re.finditer(ss_seq, seq[seq_s:seq_e].upper())])
            
                for ss_3p in alt_3p_ss:
                    if ss_3p in annotated_ss_3p: 
                        source = "NAGNNAG_annotated"
                    else:
                        source = "NAGNNAG"
                     
                    exon_paths = {"A":[0,1], "B":[0,2]}
                    NAG_exon = strand_d == "FWD" and [ss_3p, ds_exon[1]] or [ds_exon[0], ss_3p]
                    EXONS  = [us_exon, NAG_exon, ds_exon]
                    FEATURE_ID = ",".join(["%s:%d-%d"%(self.contig, e[0], e[1]) for e in EXONS])
                    GENE_NAME = ",".join(["%s"%gi['gene_name'] for gi in gene_inf])
                    GENE_ID = ",".join(["%s"%gi['gene_ID'] for gi in gene_inf])
                    G_START = min(us_exon[0], ds_exon[0])
                    G_END = max(us_exon[1], ds_exon[1])
                    STRAND = strand_d == "FWD" and 1 or -1

                    nagNnag_T = Transcript(contig = self.contig, 
                                           feature_ID = FEATURE_ID,
                                           exons = EXONS, 
                                           gene_name = GENE_NAME,
                                           gene_ID = GENE_ID,
                                           g_start = G_START,
                                           g_end = G_END,
                                           strand = STRAND)
                     
                    gff_s = nagNnag_T.gff_string(exon_paths, source)
                    bed_s = nagNnag_T.bed_string(exon_paths, source, True)
                    F_gff.write(gff_s)
                    F_bed.write(bed_s)

    def enumerate_splice_junctions(self, seq):
        fwd_5p = defaultdict(int)
        fwd_3p = defaultdict(int)
        rev_5p = defaultdict(int)
        rev_3p = defaultdict(int)

        for F_5p, F_3p_list in self.F_5p_3p_ss.iteritems():
            for F_3p in F_3p_list: 
                fwd_3p[seq[F_3p-2:F_3p]]+=1
        for F_3p, F_5p_list in self.F_3p_5p_ss.iteritems():
            for F_5p in F_5p_list: 
                fwd_5p[seq[F_5p:F_5p+2]]+=1
        for R_5p, R_3p_list in self.R_5p_3p_ss.iteritems():
            for R_3p in R_3p_list: 
                rev_3p[revcomp(seq[R_3p:R_3p+2])]+=1
        for R_3p, R_5p_list in self.R_3p_5p_ss.iteritems():
            for R_5p in R_5p_list: 
                rev_5p[revcomp(seq[R_5p-2:R_5p])]+=1 
        
        print "fwd"
        print fwd_5p, np.sum(fwd_5p.values())
        print fwd_3p, np.sum(fwd_3p.values())
        print "rev"
        print rev_5p, np.sum(rev_5p.values())
        print rev_3p, np.sum(rev_3p.values())


def init_splice_graphs_from_gff3(fn_gff, **kwargs):
    """
    limit_info={"gff_id": ["13"]}
    """
    contigs = kwargs.get('contigs', [])
    features = kwargs.get('features', ["transcript", "mRNA", "protein_coding"])

    limit_info = {"gff_id": contigs}
    GFF_parser = GFF.GFFParser()
    
    SGs_by_contig = {}
    for rec in GFF_parser.parse_in_parts(open(fn_gff), limit_info=limit_info):
        #iterate over chromosomes
        contig = rec.id
        print >>stderr, "parsing records for record id:%s..."%rec.id
        #contig_seq = fa.get_sequence(rec.id)
                   
        F_3p_5p_ss = {} 
        F_5p_3p_ss = {}  
        R_3p_5p_ss = {} 
        R_5p_3p_ss = {} 

        F_exon_s_e = {} 
        F_exon_e_s = {} 
        R_exon_s_e = {} 
        R_exon_e_s = {} 

        F_5p_to_gene_info = {}
        F_3p_to_gene_info = {}
        R_5p_to_gene_info = {}
        R_3p_to_gene_info = {}

        for feature in rec.features:
            if feature.type in features:
                t = Transcript.init_from_feature(contig, feature)
                #alt_ss = t.get_all_3pSS(contig_seq)
                if t.strand == 1:
                    t.get_splice_junctions(F_3p_5p_ss, 
                                           F_5p_3p_ss, 
                                           F_exon_s_e, 
                                           F_exon_e_s, 
                                           F_5p_to_gene_info, 
                                           F_3p_to_gene_info)
                else:
                    t.get_splice_junctions(R_3p_5p_ss, 
                                           R_5p_3p_ss, 
                                           R_exon_s_e, 
                                           R_exon_e_s, 
                                           R_5p_to_gene_info, 
                                           R_3p_to_gene_info)
        
        SGs_by_contig[contig] = SpliceGraph(contig, F_3p_5p_ss = F_3p_5p_ss, 
                                                    F_5p_3p_ss = F_5p_3p_ss, 
                                                    R_3p_5p_ss = R_3p_5p_ss, 
                                                    R_5p_3p_ss = R_5p_3p_ss, 
                                                    F_exon_s_e = F_exon_s_e,
                                                    F_exon_e_s = F_exon_e_s,
                                                    R_exon_s_e = R_exon_s_e,
                                                    R_exon_e_s = R_exon_e_s,
                                                    F_5p_to_gene_info = F_5p_to_gene_info,
                                                    F_3p_to_gene_info = F_3p_to_gene_info,
                                                    R_5p_to_gene_info = R_5p_to_gene_info,
                                                    R_3p_to_gene_info = R_3p_to_gene_info)
    return SGs_by_contig
    
#for rec in GFF_parser.parse_in_parts(open(o.fn_gff), limit_info=limit_info):
#    #iterate over chromosomes
#    contig_seq = fa.get_sequence(rec.id)
#    
#    D_to_A, A_to_D = get_donor_acceptor_pairs(rec)
#    for feature in rec.features:
#        if feature.type in ["mRNA", "transcript"]:
#            t = transcript.init_from_feature(feature)
#            alt_ss = t.get_all_3pSS(contig_seq)
#            print alt_ss
