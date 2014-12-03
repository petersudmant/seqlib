import argparse
import pandas as pd
import pdb
from sys import stderr
import numpy as np
from junction_counter import JunctionCounter
from fastahack import FastaHack
import pysam

def get_nagnags(all_juncs):
    plus = get_nagnags_by_strand(all_juncs, "+")
    minus = get_nagnags_by_strand(all_juncs, "-")
    all_nagnags = sorted(plus+minus, key=lambda x: (x['contig'],x['ss_5p']))
    return all_nagnags


def get_nagnags_by_strand(all_juncs, strand):
    
    all_juncs = all_juncs.copy()
    all_juncs = all_juncs.drop_duplicates()
    if strand == "-": 
        ss_5p_key = "j_right" 
        ss_3p_key = "j_left"
    else:
        ss_5p_key = "j_left" 
        ss_3p_key = "j_right"

    juncs = all_juncs[all_juncs["strand"]==strand]
    juncs = juncs.sort(['contig',ss_5p_key, ss_3p_key])

    prev_ss_5p = None
    prev_ss_3p = None
    
    nagnags = []

    for i, row in juncs.iterrows():
        ss_5p = row[ss_5p_key]
        ss_3p = row[ss_3p_key]
        contig = row['contig']
        if prev_ss_5p and ss_5p == prev_ss_5p:
            #print contig, ss_5p, ss_3p, prev_ss_3p, 
            if strand == "+":
                j_left1, j_right1 = ss_5p, prev_ss_3p
                j_left2, j_right2 = ss_5p, ss_3p
                ss_3p_prox, ss_3p_dist = prev_ss_3p, ss_3p
            else:
                j_left1, j_right1 = prev_ss_3p, ss_5p
                j_left2, j_right2 = ss_3p, ss_5p 
                ss_3p_prox, ss_3p_dist = ss_3p, prev_ss_3p 

            size = abs(ss_3p-prev_ss_3p)
            
            nagnags.append({"contig":contig,
                            "ss_5p":ss_5p,
                            "ss_3p_prox": ss_3p_prox,
                            "ss_3p_dist": ss_3p_dist,
                            "strand": strand})

        else: 
            prev_ss_5p = ss_5p
            prev_ss_3p = ss_3p
    
    return nagnags

        
def get_psi(nagnag, bam, readlen, min_overhang):
    
    contig = nagnag['contig']
    if nagnag['strand'] == "+":
        j_prox = tuple([nagnag['ss_5p'], nagnag['ss_3p_prox']])
        j_dist = tuple([nagnag['ss_5p'], nagnag['ss_3p_dist']])
        s = nagnag['ss_5p']
        e = nagnag['ss_3p_dist']
    else:
        j_prox =  tuple([nagnag['ss_3p_prox'], nagnag['ss_5p']])
        j_dist =  tuple([nagnag['ss_3p_dist'], nagnag['ss_5p']])
        s = nagnag['ss_3p_dist']
        e = nagnag['ss_5p']
    
    j_prox_w_contig = tuple([contig, j_prox[0], j_prox[1]])
    j_dist_w_contig = tuple([contig, j_dist[0], j_dist[1]])

    junction_inf = {j_prox: JunctionCounter(j_prox_w_contig, readlen),
                    j_dist: JunctionCounter(j_dist_w_contig, readlen)}

    #print junction_inf.keys()
    for read in bam.fetch(reference=contig, start = s, end = e):
        curr_read_pos = 0
        for i in xrange(1,len(read.blocks)):
            read_junc = tuple([read.blocks[i-1][1]-1, read.blocks[i][0]])
            curr_read_pos += read.blocks[i-1][1]-read.blocks[i-1][0]
            if read_junc in junction_inf: 
                junction_inf[read_junc].add_read(curr_read_pos)
        del read

    prox_n = junction_inf[j_prox].n_supporting_reads(min_overhang) 
    dist_n = junction_inf[j_dist].n_supporting_reads(min_overhang) 
    
    if dist_n+prox_n == 0:
        psi = 0
    else:
        psi = prox_n/(dist_n+prox_n)

    return psi, prox_n, dist_n

def get_readlen(bam):
    for read in bam.fetch():
        return read.rlen

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_fasta")
    parser.add_argument("--fn_juncs")
    parser.add_argument("--fn_bam")
    parser.add_argument("--fn_output")
    parser.add_argument("--min_overhang", type=int, default=8)
    o = parser.parse_args()
    
    FOUT = open(o.fn_output, 'w')

    fa = FastaHack(o.fn_fasta)
    bam = pysam.Samfile(o.fn_bam, 'rb')
    contigs = [contig for contig in fa.names]
    #contigs = ["chr4"]
    #contigs = ['chrUn.004.790']

    junc_table = pd.read_csv(o.fn_juncs, 
                             header=None, 
                             delimiter="\t",
                             names=["contig",
                                    "j_left",
                                    "j_right",
                                    "strand"])
    nagnags = get_nagnags(junc_table)
    """
    Vectors are readlen -1 because you need at least one base overlapping the
    junction
    """
    
    readlen = get_readlen(bam)
    total_read_count_vect = np.zeros(readlen-1)
    FOUT.write('contig\tss_5p\tss_3p_prox\tss_3p_dist\tsize\tpsi\tprox_count\tdist_count\n') 
    for i, nagnag in enumerate(nagnags):
        if i %1000 == 0:
            print i
        psi, prox_n, dist_n = get_psi(nagnag, bam, readlen, o.min_overhang)
        nagnag['psi'] = psi
        size = abs(nagnag['ss_3p_prox']-nagnag['ss_3p_dist'])
        FOUT.write('{contig}\t{ss_5p}\t{ss_3p_prox}\t{ss_3p_dist}\t{size}\t{psi}'
                                                   '\t{prox_n}\t{dist_n}\n'.format(contig = nagnag['contig'],
                                                                                   psi = psi,
                                                                                   size=size,
                                                                                   prox_n = prox_n,
                                                                                   dist_n = dist_n,
                                                                                   ss_5p =nagnag['ss_5p'],
                                                                                   ss_3p_prox = nagnag['ss_3p_prox'],
                                                                                   ss_3p_dist = nagnag['ss_3p_dist']))

