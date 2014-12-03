import argparse
import pdb
import time
import numpy as np
from fastahack import FastaHack
from collections import defaultdict
import BCBio.GFF as GFF
import collections

from bx.intervals.cluster import ClusterTree
import time
from junction_writer import JunctionWriter
import pysam

from guppy import hpy
from junction_counter import JunctionCounter

#was calculated using the following equations:
#    pi = reads at offset i / total reads to junction window
#    Entopy = - sumi(pi * log(pi) / log2)

def get_junction_entropy(contig, start, end, juncs, bam, readlen, entropy_min_overhang, FOUT):
    
    ##!!! NEED TO DO SOMETHING WHEN NOTHIN IN THE FILE!
    
    #for pcol in bam.pileup(reference=contig, start = start, end = end):
    read_count_vect = np.zeros(readlen-1)
    junction_inf = {}
    strand = juncs[0][-1]

    for junc in juncs:
        junction_inf[tuple([junc[1],junc[2]])] = JunctionCounter(junc, readlen)
    
    #print junction_inf.keys()
    for read in bam.fetch(reference=contig, start = start, end = end):
        read_juncs = []
        curr_read_pos = 0
        for i in xrange(1,len(read.blocks)):
            read_junc = tuple([read.blocks[i-1][1]-1, read.blocks[i][0]])
            curr_read_pos += read.blocks[i-1][1]-read.blocks[i-1][0]
            if read_junc in junction_inf: 
                junction_inf[read_junc].add_read(curr_read_pos)
        del read

    for junc, inf in junction_inf.iteritems():
        e = inf.calc_entropy(entropy_min_overhang)
        read_count_vect+=inf.read_counts
        min_overhang = inf.min_lr_overhang()
        FOUT.write("{contig}\t{j_l}\t{j_r}\t{strand}\t"
                    "{entropy}\t{min_overhang}\n".format(contig=contig,
                                                         j_l = junc[0],
                                                         j_r = junc[1],
                                                         strand = strand,
                                                         entropy=e,
                                                         min_overhang = min_overhang))
        del inf
        
    return read_count_vect
#import re
#import operator
def get_juncs(tbx, contig):
    if not contig in tbx.contigs:
        return

    for c in tbx.fetch(contig, parser=pysam.asTuple()):
        yield c[0], int(c[1]), int(c[2]), c[3]

def get_cluster_trees(tbx_juncs, contigs):
    
    junc_clusts_by_strand_contigs = {"+":collections.defaultdict(lambda:
                                        ClusterTree(0,2)), 
                                     "-" :collections.defaultdict(lambda:
                                        ClusterTree(0,2))}
    id=0
    juncs_by_id = {}
    for contig in contigs:
        for junc in get_juncs(tbx_juncs, contig):
            contig, s, e, strand = junc
            junc_clusts_by_strand_contigs[strand][contig].insert(s,e,id)
            juncs_by_id[id] = junc
            id+=1

    return dict(junc_clusts_by_strand_contigs), juncs_by_id

def get_readlen(bam):
    for read in bam.fetch():
        return read.rlen

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_fasta")
    parser.add_argument("--fn_juncs")
    parser.add_argument("--fn_bam")
    parser.add_argument("--fn_output")
    parser.add_argument("--entropy_min_overhang", type=int, default=8)
    o = parser.parse_args()
    
    FOUT = open(o.fn_output, 'w')
    FOUT_pileup = open("%s.pileup"%o.fn_output, 'w')

    fa = FastaHack(o.fn_fasta)
    tbx_juncs = pysam.Tabixfile(o.fn_juncs) 
    bam = pysam.Samfile(o.fn_bam, 'rb')
    contigs = [contig for contig in fa.names]
    #contigs = ["chr4"]
    #contigs = ['chrUn.004.790']
    junc_cluster_trees, juncs_by_id = get_cluster_trees(tbx_juncs, contigs)
    """
    Vectors are readlen -1 because you need at least one base overlapping the
    junction
    """
    
    readlen = get_readlen(bam)
    total_read_count_vect = np.zeros(readlen-1)
    junction_clusts_assessed = 0 
    for strand, junc_clusters_by_contig in junc_cluster_trees.iteritems():
        print strand
        for contig, junc_clusters in junc_clusters_by_contig.iteritems():
            print contig
            t=time.time()
            for start, end, junction_ids in junc_clusters.getregions():
                junction_clusts_assessed += 1 
                #print start, end, junction_ids, [juncs_by_id[i] for i in junction_ids]
                juncs = [juncs_by_id[i] for i in junction_ids]
                total_read_count_vect += get_junction_entropy(contig, start, end, juncs, bam, readlen, o.entropy_min_overhang, FOUT)
            #print "time for {strand} {contig} {t}".format(contig=contig, 
            #                                             strand = strand,
            #                                             t=time.time()-t)
            #print total_read_count_vect
            #if junction_clusts_assessed>100: 
            #    h = hpy()
            #    print h.heap()
            #    junction_clusts_assessed = 0
            #    pdb.set_trace()
    FOUT_pileup.write("position\tread_count\n")
    for i in xrange(len(total_read_count_vect)):
        FOUT_pileup.write("{pos}\t{count}\n".format(pos=i,count=total_read_count_vect[i]))
        

