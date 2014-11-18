import argparse
import gzip

if __name__=="__main__":


    parser = argparse.ArgumentParser()
    parser.add_argument("--fastq1")
    parser.add_argument("--fastq2")
    parser.add_argument("--n_reads",default=1000000,type=int)
    parser.add_argument("--gzipped",default=False, action='store_true')
    o = parser.parse_args() 

    if o.gzipped:
        F1 = gzip.open(o.fastq1)
        F2 = gzip.open(o.fastq2)
        for i in xrange(o.n_reads):
            l1_1 = F1.readline()
            l1_2 = F1.readline()
            l1_3 = F1.readline()
            l1_4 = F1.readline()
            l2_1 = F2.readline()
            l2_2 = F2.readline()
            l2_3 = F2.readline()
            l2_4 = F2.readline()
            print("%s%s%s%s"%(l1_1,
                              l1_2,
                              l1_3,
                              l1_4.rstrip()))
            print("%s%s%s%s"%(l2_1,
                              l2_2,
                              l2_3,
                              l2_4.rstrip()))
        
