import argparse
import pandas as pd
import pdb
from sys import stderr
import numpy as np
from JunctionCounter import JunctionCounter

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

    junction_inf = {j_prox: JunctionCounter(j_prox, readlen),
                    j_dist: JunctionCounter(j_dist, readlen)}

    #print junction_inf.keys()
    for read in bam.fetch(reference=contig, start = s, end = e):
        if len(read.blocks)==2:
            read_junc = tuple([read.blocks[0][1]-1, read.blocks[1][0]])
            if read_junc in junction_inf: 
                junction_inf[read_junc].add_read(read.pos)
        del read
    
    prox_n = junction_inf[j_prox].n_supporting_reads(min_overhang) 
    dist_n = junction_inf[j_dist].n_supporting_reads(min_overhang) 
    
    return prox_n/(dist_n+prox_n)

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
    tbx_juncs = pysam.Tabixfile(o.fn_juncs) 
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
                                    "strand")
    nagnags = get_nagnags(junc_table)
    """
    Vectors are readlen -1 because you need at least one base overlapping the
    junction
    """
    
    readlen = get_readlen(bam)
    total_read_count_vect = np.zeros(readlen-1)
    FOUT.write('contig\tss_5p\tss_3p_prox\tss_3p_dist\tpsi\n') 
    for nagnag in nagnags:
        psi = get_psi(nagnag, bam, readlen, o.min_overhang)
        nagnag['psi'] = psi
        FOUT.write('{contig}\t{ss_5p}\t{ss_3p_prox}\t{ss_3p_dist}\t{psi}\n'.format(nagnag)) 

