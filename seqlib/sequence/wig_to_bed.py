import pdb
#from  bx.arrays.wiggle import WiggleReader
#from bx import wiggle
from bx.arrays import wiggle
import gzip
from fastahack import FastaHack
import numpy as np
import argparse
#help(wig)


class np_wiggle(object):
    def __init__(self, path, contig_len, gzipped=True):
        assert gzipped==True
        self.values = np.zeros(contig_len)
        w  = wiggle.WiggleReader(gzip.GzipFile(path))
        #w  = wiggle.WiggleReader(open(path))
        counter = 0
        for l in w:
            counter+=1
            print "\t".join([str(i) for i in l])

if __name__=="__main__":
    parser = argparse()
    parser.add_argument("--fn_wig")
    fasta = "/net/uorf/data/psudmant/genomes/UCSC/UCSC/mus_musculus/mm9/fasta/mm9.fa" 
    phastConsDir = "/net/uorf/data/psudmant/genomes/UCSC/UCSC/mus_musculus/mm9/phastCons30way/vertebrate/"
    
    fa = FastaHack(fasta)
    for contig in fa.names:
        #if contig != "chr1": continue
        l = fa.get_sequence_length(contig)
        print contig
        w = np_wiggle("%s/%s.data.gz"%(phastConsDir,contig), l)
        #w = np_wiggle("%s/%s.data"%(phastConsDir,contig), l)
        pdb.set_trace()



