import argparse
import pdb
import os
import numpy as np
import pandas as pd
from collections import defaultdict


def m_score(deltas):

    s = np.sum(deltas>0)-np.sum(deltas<0)

    t = np.tile(deltas,(1000,1))
    rsign = ((np.random.rand(1000,deltas.shape[0])>0.5)*2)-1
    rnd_deltas = t*rsign
    dist = np.sum(rnd_deltas>0,1)-np.sum(rnd_deltas<0,1)
    mu = np.mean(dist)
    sd = np.std(dist)
    
    MZ=(s-mu)/sd
    return s, MZ

def get_M_scores(sample_order, t_by_comparison, events, BF_min):

    comparisons = t_by_comparison.keys()
    dvect = [sample_order.index(c[0])<sample_order.index(c[1]) for c in comparisons]
    assert np.all(np.array(dvect))
    
    m_table = []
    
    for i,e in enumerate(events):
        BFs = []
        deltas = []
        for c in comparisons:
            if t_by_comparison[c].loc[e]['bayes_factor']>=BF_min:
                deltas.append(t_by_comparison[c].loc[e]['sample2_posterior_mean']-t_by_comparison[c].loc[e]['sample1_posterior_mean'])
        BFs = np.array(BFs)
        deltas = np.array(deltas)

        if deltas.shape[0]>1:
            m, MZ = m_score(deltas)
        else:
            m, MZ = 0, 0
        
        m_table.append({"compressed_ID":e, "M":m, "MZ":MZ})
    
    t = pd.DataFrame(m_table)
    return t


if __name__=="__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("--BF_cutoff", type=float, default=5)
    parser.add_argument("--fn_MISO_comparisons", nargs="+")
    parser.add_argument("--sample_order", nargs="+")
    parser.add_argument("--MISO_ids")
    parser.add_argument("--fn_out")
    o = parser.parse_args()
    
    T_ids = pd.read_csv(o.MISO_ids, header=0,sep="\t")
    T_ids.set_index("compressed_ID")

    t_by_comparison = {}
    events = None
    
    for fn_compare in o.fn_MISO_comparisons:
        comp = fn_compare.split("/")[-2]
        s1, s2 = comp.split("_vs_")
        t=pd.read_csv(fn_compare, header=0, sep="\t")
        t = t[['event_name','bayes_factor','sample1_posterior_mean','sample2_posterior_mean']]
        """
        remove any dudes with ","s, ie, non-binary events
        """
        t = t[t.apply(lambda x: not "," in str(x['bayes_factor']), axis=1)]
        t['sample1_posterior_mean'] = t['sample1_posterior_mean'].astype("float")
        t['sample2_posterior_mean'] = t['sample2_posterior_mean'].astype("float")
        t['bayes_factor'] = t['bayes_factor'].astype("float")
        t=t.set_index(['event_name'])
        t_by_comparison[tuple([s1,s2])] = t
        if events:
            events = events.intersection(set(t.index))
        else:
            events = set(t.index)
    
    psi_T = pd.DataFrame()
    for c, t in t_by_comparison.items():
        s1, s2 = c
        if not psi_T.empty:
            if not s1 in psi_T.columns:
                j = t['sample1_posterior_mean']
                j.name = s1
                psi_T = psi_T.join(j)
            if not s2 in psi_T.columns:
                j = t['sample2_posterior_mean']
                j.name = s2
                psi_T = psi_T.join(j)
        else:
            psi_T=t[['sample1_posterior_mean','sample2_posterior_mean']]
            psi_T.columns = [s1,s2]
        
    T_MZ = get_M_scores(o.sample_order, t_by_comparison, events, o.BF_cutoff)
    
    """now, join T_MZ"""

    T_MZ = pd.merge(T_MZ, T_ids, left_on="compressed_ID", right_on="compressed_ID")
    T_MZ = pd.merge(T_MZ, psi_T, left_on="compressed_ID", right_index=True)
    header=["compressed_ID","ID", "type", "M", "MZ"]+o.sample_order
    T_MZ.to_csv(o.fn_out, columns=header,index=False, sep="\t")





