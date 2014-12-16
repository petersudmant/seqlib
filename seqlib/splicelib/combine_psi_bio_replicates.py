import argparse
import pdb
import time
import numpy as np
from fastahack import FastaHack
from collections import defaultdict
import BCBio.GFF as GFF
from junction_writer import JunctionWriter
import pysam
import pandas as pd
from sys import stderr
from scipy.stats import pearsonr

def assess_bio_replicates(T, fn_out):

    tissues = list(np.unique(T["tissue"].values))
    
    grouped = T.groupby(["contig","ss_5p","ss_3p_prox","ss_3p_dist","size", "tissue"])
    nn_count_summary = grouped['n_reads'].agg({"n_reads_mean":np.mean,"n_reads_min":np.min})
    T = T.join(nn_count_summary, on=["contig","ss_5p","ss_3p_prox","ss_3p_dist","size","tissue"])
    T = T[T['n_reads_min']>=10]
    
    bio_rep_sum = []
    for tissue in tissues:
        t_T = T[T['tissue']==tissue]
        replicates = list(np.unique(t_T["sample"].values))
        
        for i in range(len(replicates)):
            for j in range(i, len(replicates)):
                if i==j: continue
                rep1, rep2 = replicates[i], replicates[j]
                #c = np.corrcoef(t_T[t_T["sample"]==rep1]["psi"],t_T[t_T["sample"]==rep2]["psi"])
                r,p = pearsonr(t_T[t_T["sample"]==rep1]["psi"],t_T[t_T["sample"]==rep2]["psi"])
                bio_rep_sum.append({"r.squared":r*r,"tissue":tissue, "sample1":rep1, "sample2":rep2})

    
    df = pd.DataFrame(bio_rep_sum)
    df.to_csv(fn_out, sep="\t", index=False)

def get_sum_table(T, sample_tissue_tups, fn_out):
    
    #nn_summary = T.groupby(["contig","ss_5p","ss_3p_prox","ss_3p_dist","size", "tissue"]).agg([np.mean,np.std])
    #nn_summary = grouped.agg({"psi":[np.mean,np.std],'prox_count':np.mean,"dist_count":np.mean}, as_index=False)
    #nn_nreads_summary = grouped['psi'].agg({"psi_mean":np.mean,"psi_std":np.std})
    grouped = T.groupby(["contig","ss_5p","ss_3p_prox","ss_3p_dist","size", "tissue"])
    nn_psi_summary = grouped['psi'].agg({"psi_mean":np.mean,"psi_std":np.std})
    nn_count_summary = grouped['n_reads'].agg({"n_reads_mean":np.mean,"n_reads_min":np.min})
    nn_psi_summary = nn_psi_summary.join(nn_count_summary).reset_index()
    
    switch_scores = nn_psi_summary.groupby(["contig","ss_5p","ss_3p_prox","ss_3p_dist","size"])['psi_mean'].agg({"switch_score":lambda x:np.amax(x)-np.amin(x)}) 
    nn_psi_summary = nn_psi_summary.join(switch_scores, on=["contig","ss_5p","ss_3p_prox","ss_3p_dist","size"])

    nn_psi_summary.to_csv(fn_out,sep="\t", index=False)
    #for tissue in tissues:
        
def read_tables(fns, l_get_sample_tissue):
                
    
    tables_by_sample_tissue = {}
    sample_tissue_tups = [] 
    for f in fns:
        sample, tissue = l_get_sample_tissue(f) 
        sample_tissue = "%s_%s"%(sample,tissue)
        sample_tissue_tup = tuple([sample,tissue])
        sample_tissue_tups.append(sample_tissue_tup)
        
        t = pd.read_csv(f, 
                        header=0, 
                        delimiter="\t")
                        #names = ["contig",
                        #         "ss_5p",
                        #         "ss_3p_prox",
                        #         "ss_3p_dist",
                        #         "size",
                        #         "psi",
                        #         "prox_count",
                        #         "dist_count",
                        #         ])
        t['tissue'] = tissue
        t['sample'] = sample
        tables_by_sample_tissue[sample_tissue_tup] = t
        print >> stderr, tissue   
    
    tables = tables_by_sample_tissue.values() 
    T = pd.concat(tables)
    return T, sample_tissue_tups

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_out_all")
    parser.add_argument("--fn_out_summary")
    parser.add_argument("--fn_out_bio_replicate_summary")
    parser.add_argument("--sample_tissue_lambda", default="""lambda x: x.split("/")[-1].split(".")[0].split("_")""")
    parser.add_argument("--fn_inputs", nargs="*")
    o = parser.parse_args()
    
    sample_tissue_L = eval(o.sample_tissue_lambda) 
    T, sample_tissue_tups = read_tables(o.fn_inputs, sample_tissue_L)
    T["n_reads"] = T["prox_count"]+T["dist_count"]
    T.to_csv(o.fn_out_all,index=False, sep="\t")
    assess_bio_replicates(T, o.fn_out_bio_replicate_summary)
    sum_table = get_sum_table(T, sample_tissue_tups, o.fn_out_summary)



