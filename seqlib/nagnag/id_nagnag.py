import BCBio.GFF as GFF
import argparse
import pdb
import time
import numpy as np
from fastahack import FastaHack
import re
import string 
from collections import defaultdict
import operator

trans = string.maketrans('ATCGatcg', 'TAGCtacg')
def revcomp(s):
    return s.translate(trans)[::-1]

class transcript:

    def __init__(self, **kwargs):
        self.feature_ID = kwargs['feature_ID'] 
        self.gene_name = kwargs['gene_name'] 
        self.gene_ID = kwargs['gene_ID'] 
        self.g_start = kwargs['g_start'] 
        self.g_end = kwargs['g_end'] 
        self.exons = kwargs['exons']
        self.strand = kwargs['strand']
    
    def print_out(self):
        print self.feature_ID
        print "\t", self.gene_name
        print "\t", self.gene_ID 
        print "\t", self.g_start
        print "\t", self.g_end
        print "\t", self.exons
        print "\t", self.strand
    
    def get_all_3pSS(self, seq):

        for i in xrange(len(self.exons)-1):
            e1_e = self.exons[i][1]
            e2_s = self.exons[i+1][0]
            
            if self.strand == 1:
                alts = np.array([m.start() for m in re.finditer("AG",seq[e1_e:e2_s])])
            else:
                alts = np.array([m.start() for m in re.finditer("CT",seq[e1_e:e2_s])])
            return alts 

    def get_splice_junctions(self, seq):
        ret = []
        for i in xrange(len(self.exons)-1):
            e1_e = self.exons[i][1]
            e2_s = self.exons[i+1][0]

            if self.strand == 1:
                ret.append(tuple([seq[e1_e:e1_e+2],seq[e2_s-2:e2_s]]))
            else:
                ret.append(tuple([revcomp(seq[e2_s-2:e2_s]), revcomp(seq[e1_e:e1_e+2])])) 
            
        return ret

    @classmethod
    def init_from_feature(cls, feature):
        kwargs = {
                  'feature_ID': feature.id,
                  'gene_name': feature.qualifiers['gene_name'][0],
                  'gene_ID': feature.qualifiers['geneID'][0],
                  'strand': feature.strand,
                  'g_start': feature.location.start.position,
                  'g_end': feature.location.end.position,
                  'exons':[[s.location.start.position,
                            s.location.end.position] for s in feature.sub_features]
                  }
                  
        return cls(**kwargs)

def get_donor_acceptor_pairs(rec):
    D_to_A, A_to_D = 

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_fasta")
    parser.add_argument("--fn_gff")
    o = parser.parse_args()
    
    fa = FastaHack(o.fn_fasta)

    GFF_parser = GFF.GFFParser()
    limit_info={"gff_id": ["13"]}
    #limit_info={}
    
    for rec in GFF_parser.parse_in_parts(open(o.fn_gff), limit_info=limit_info):
        #iterate over chromosomes
        contig_seq = fa.get_sequence(rec.id)
        
        D_to_A, A_to_D = get_donor_acceptor_pairs(rec)
        for feature in rec.features:
            if feature.type in ["mRNA", "transcript"]:
                t = transcript.init_from_feature(feature)
                alt_ss = t.get_all_3pSS(contig_seq)
                print alt_ss
