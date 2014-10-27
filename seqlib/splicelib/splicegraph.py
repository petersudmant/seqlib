import numpy as np
import pdb
import argparse
import BCBio.GFF as GFF
from transcript import Transcript
from collections import defaultdict
from sys import stderr

import string 
trans = string.maketrans('ATCGatcg', 'TAGCtacg')

def revcomp(s):
    return s.translate(trans)[::-1]

class SpliceGraph(object):

    def __init__(self, contig, F_3p_5p_ss, F_5p_3p_ss, R_3p_5p_ss, R_5p_3p_ss):
        
        """
          5'      3'
        + GT------AG
        
          3'      5'
        - CT------AC
        """
        
        self.contig = contig

        self.F_3p_5p_ss = F_3p_5p_ss 
        self.F_5p_3p_ss = F_5p_3p_ss 
        
        self.R_3p_5p_ss = R_3p_5p_ss
        self.R_5p_3p_ss = R_5p_3p_ss
    
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

def init_splice_graphs_from_gff(fn_gff, contigs = [], features = ["transcript", "mRNA"]):
    """
    limit_info={"gff_id": ["13"]}
    """
    limit_info = {"gff_id": contigs}
    GFF_parser = GFF.GFFParser()
    
    SGs_by_contig = {}
    
    for rec in GFF_parser.parse_in_parts(open(fn_gff), limit_info=limit_info):
        #iterate over chromosomes
        contig = rec.id
        print >>stderr, "parsing records for record id:%s..."%rec.id
        #contig_seq = fa.get_sequence(rec.id)
                   
        F_3p_5p_ss = defaultdict(list) 
        F_5p_3p_ss = defaultdict(list)
        R_3p_5p_ss = defaultdict(list)
        R_5p_3p_ss = defaultdict(list)

        F_exon_5p_3p = defaultdict(list) 
        F_exon_3p_5p = defaultdict(list)
        R_exon_5p_3p = defaultdict(list)
        R_exon_3p_5p = defaultdict(list)

        for feature in rec.features:
            if feature.type in features:
                t = Transcript.init_from_feature(feature)
                #alt_ss = t.get_all_3pSS(contig_seq)
                if t.strand == 1:
                    t.get_splice_junctions(F_3p_5p_ss, F_5p_3p_ss, F_exon_3p_5p, F_exon_5p_3p)
                else:
                    t.get_splice_junctions(R_3p_5p_ss, R_5p_3p_ss, R_exon_3p_5p, R_exon_5p_3p)
        
        SGs_by_contig[contig] = SpliceGraph(contig, F_3p_5p_ss, F_5p_3p_ss, R_3p_5p_ss, R_5p_3p_ss)
                
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
