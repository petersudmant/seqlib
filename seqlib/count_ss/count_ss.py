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
from sys import stderr

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
    
    def print_splice_junctions(self, seq):

        for i in xrange(len(self.exons)-1):
            e1_e = self.exons[i][1]
            e2_s = self.exons[i+1][0]
            #print seq[e1_e:e1_e+2],seq[e2_s-2:e2_s]
            #print seq[e1_e:e1_e+2],seq[e2_s-2:e2_s]

            if self.strand == 1:
                print seq[e1_e:e1_e+2],seq[e2_s-2:e2_s]
                print seq[e1_e:e1_e+2],seq[e2_s-2:e2_s]
                alts = np.array([m.start() for m in re.finditer("AG",seq[e1_e:e2_s])])
            else:
                print revcomp(seq[e2_s-2:e2_s]), revcomp(seq[e1_e:e1_e+2]) 
                print revcomp(seq[e2_s-2:e2_s]), revcomp(seq[e1_e:e1_e+2]) 
                alts = np.array([m.start() for m in re.finditer("CT",seq[e1_e:e2_s])])
            #print alts

    def get_splice_junctions(self, seq, count_juncs):
        ret = []
        for i in xrange(len(self.exons)-1):
            e1_e = self.exons[i][1]
            e2_s = self.exons[i+1][0]

            if self.strand == 1:
                ret.append(tuple([seq[e1_e:e1_e+2],seq[e2_s-2:e2_s]]))
            else:
                ret.append(tuple([revcomp(seq[e2_s-2:e2_s]), revcomp(seq[e1_e:e1_e+2])])) 
            count_juncs[ret[-1]]+=1
            
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


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_fasta")
    parser.add_argument("--fn_gff")
    o = parser.parse_args()
    
    fa = FastaHack(o.fn_fasta)

    GFF_parser = GFF.GFFParser()
    limit_info={"gff_id": [str(i) for i in xrange(1,23)]+["X","Y"]}
    #limit_info={"gff_id": ["20"]}
    
    count_juncs = defaultdict(int)

    for rec in GFF_parser.parse_in_parts(open(o.fn_gff), limit_info=limit_info):
        contig_seq = fa.get_sequence(rec.id)
        print >>stderr, rec.id
        for feature in rec.features:
            if feature.type in ["mRNA", "transcript"]:
                t = transcript.init_from_feature(feature)
                t.get_splice_junctions(contig_seq, count_juncs)
    
    print "sequence\tfrequency"
    for s in sorted(count_juncs.items(), key=operator.itemgetter(1)):
        print "\t".join(["%s-%s"%(s[0][0], s[0][1]),  str(s[1])])
    
