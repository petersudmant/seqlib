import tables
import numpy as np
import pandas as pd
import argparse
import sys
import pdb
import scipy.stats as stats

from coveragedata import CoverageData 


class h5FullGeneCvg_writer(object):
    
    def __init__(self, fn):
        self.h5 = tables.openFile(fn, mode='w')
        
        filt = tables.Filters(complevel=5, complib='blosc')
        self.a_cvg = self.h5.createEArray(self.h5.root, 
                                          'cvg', 
                                          tables.Float32Atom(), 
                                          shape=(0,), 
                                          filters=filt,
                                          expectedrows=65000000)
        
        self.a_idx = self.h5.createEArray(self.h5.root, 
                                          'idx', 
                                          tables.UInt16Atom(), 
                                          shape=(0,), 
                                          filters=filt,
                                          expectedrows=65000000)
        
        self.a_GeneID = self.h5.createEArray(self.h5.root, 
                                             'geneID', 
                                             tables.StringAtom(20),
                                             shape=(0,), 
                                             filters=filt,
                                             expectedrows=20000)
        
        self.a_start = self.h5.createEArray(self.h5.root, 
                                            'start', 
                                            tables.UInt16Atom(),
                                            shape=(0,), 
                                            filters=filt,
                                            expectedrows=20000)
        
        self.a_stop = self.h5.createEArray(self.h5.root, 
                                           'stop', 
                                           tables.UInt16Atom(),
                                           shape=(0,), 
                                           filters=filt,
                                           expectedrows=20000)
        
        self.a_length = self.h5.createEArray(self.h5.root, 
                                             'length', 
                                             tables.UInt16Atom(),
                                             shape=(0,), 
                                             filters=filt,
                                             expectedrows=20000)
        self.curr_idx = 0
        self.GeneIDs = []
        self.lengths = []
        self.starts = []
        self.stops = []
    
    def extend(self, cvg_obj):

        self.a_cvg.append(cvg_obj.RNA_cvg_view)
        self.a_idx.append(np.ones(cvg_obj.RNA_cvg_view.shape[0])*self.curr_idx)
        
        #self.GeneIDs.append(cvg_obj.g.gene_id)
        self.GeneIDs.append(cvg_obj.gene_id)
        self.lengths.append(cvg_obj.RNA_cvg_view.shape[0])
        self.starts.append(cvg_obj.RNA_cvg_view_START)
        self.stops.append(cvg_obj.RNA_cvg_view_STOP)

        self.curr_idx+=1
    
    def close(self):
        
        self.a_GeneID.append(np.array(self.GeneIDs))
        self.a_length.append(np.array(self.lengths))
        self.a_start.append(np.array(self.starts))
        self.a_stop.append(np.array(self.stops))

        self.h5.close()

class h5FullGeneCvg(object):
    
    def __init__(self, fn, fn_annot, annot_key):
        
        """
        idx = 1 idx per transcript / gene
        """
 
        sys.stderr.write("loading {fn}...".format(fn=fn))
        self.h5 = tables.openFile(fn, mode='r')
        
        self.cvg = self.h5.root.cvg[:]
        self.idx = self.h5.root.idx[:]

        self.GeneID = self.h5.root.geneID[:]
        self.length = self.h5.root.length[:]
        self.start = self.h5.root.start[:]
        self.stop = self.h5.root.stop[:]
        
        self.inf = pd.DataFrame({"GeneID":self.GeneID,
                                 "length":self.length,
                                 "start":self.start,
                                 "stop":self.stop, 
                                 "idx":np.unique(self.idx)})
        
        self.inf = self.inf.set_index("idx", drop=False)
        
        """
        if "filtering" then filter everything here...
        """
   
        self.table = pd.DataFrame({"cvg":self.cvg, "idx":self.idx})
        
        self.grouped_by_gene = self.table.groupby('idx') 
        self.full_stats = self.grouped_by_gene.agg([np.sum, np.mean, np.std, np.median]).reset_index(level=2, drop=True)
        self.full_stats.columns = ["sum", "mean", "std", "median"]
        self.full_stats['idx'] = self.full_stats.index
        
        self.annot = False
        
        if fn_annot:
            self.annot = True
            t = pd.read_csv(fn_annot, header=0, sep="\t")
            self.t_annot = pd.merge(self.inf, t, left_on = "GeneID", right_on="ENSEMBL_ID")[['idx',annot_key]]
            self.t_annot.columns = ['idx', 'annot']
        
        sys.stderr.write("done\n")
    
    def add_annotation(self, S, T):
        S['annot'] = "all"
        T = pd.merge(T, self.t_annot, left_on="idx", right_on="idx").reset_index()        
        S2 = T.groupby(['pos','annot']).describe().unstack()
        col_names = [".".join(x).rstrip(".")  for x in S2.columns] 
        S2.columns = col_names
        S2 = S2.reset_index()
        S = pd.concat([S,S2])
        return S 
 
    def get_stop_bp_meta(self, n_bp=300):
        sys.stderr.write("getting stop codon centered...")
        
        idxs, cvg_by_bp = [], [] 
        n=0
        for idx, g in self.grouped_by_gene:
            stop = self.inf.ix[idx]['stop'] 
            length = self.inf.ix[idx]['length']
         
            if stop >=n_bp and length >= stop+n_bp+1:
                cvg_by_bp.append(g.cvg[(stop-n_bp):(stop+n_bp+1)])
                idxs.append(np.ones(2*n_bp+1)*idx)
                n+=1
            
        sys.stderr.write("finished individual gene parsing...")
        
        T = pd.DataFrame({"idx":np.concatenate(idxs), 
                          "bp_cvg":np.concatenate(cvg_by_bp),
                          "pos":np.tile(np.arange(-n_bp,n_bp+1),n)})
        
        T = pd.merge(T, self.full_stats, left_on="idx", right_on="idx")
        T = T[T['mean']>0]
        
        bp_sum = pd.DataFrame(T.groupby("idx").apply(lambda x: np.sum(x['bp_cvg'])))
        bp_sum['idx'] = bp_sum.index
        bp_sum.columns = ['bp_sum', "idx"]
        T = pd.merge(T, bp_sum, left_on="idx", right_on="idx")
        T['bp_normalized'] = T['bp_cvg']/T['bp_sum']
        
        T = T[T['bp_sum']>0]
        
        S = T.groupby('pos').describe().unstack()
        col_names = [".".join(x) for x in S.columns]
        
        S.columns = col_names
        S['pos'] = S.index 
        
        if self.annot:
            S = self.add_annotation(S, T)
        
        sys.stderr.write("done\n")
        return S 
    
    def get_3p_bp_meta(self, n_bp=1000):
        sys.stderr.write("getting 3' bp meta...")
         
        idxs, cvg_by_bp = [], [] 
        n=0
        for idx, g in self.grouped_by_gene:
            if g.cvg.shape[0] >=n_bp:
                cvg_by_bp.append(g.cvg[-n_bp:])
                idxs.append(np.ones(n_bp)*idx)
                n+=1
        
        sys.stderr.write("finished individual gene parsing...")
        
        T = pd.DataFrame({"idx":np.concatenate(idxs), 
                          "bp_cvg":np.concatenate(cvg_by_bp),
                          "pos":np.tile(np.arange(-n_bp,0)+1,n)})
        
        T = pd.merge(T, self.full_stats, left_on="idx", right_on="idx")
        T = T[T['mean']>0]
        
        bp_sum = pd.DataFrame(T.groupby("idx").apply(lambda x: np.sum(x['bp_cvg'])))
        bp_sum['idx'] = bp_sum.index
        bp_sum.columns = ['bp_sum', "idx"]
        T = pd.merge(T, bp_sum, left_on="idx", right_on="idx")
        T['bp_normalized'] = T['bp_cvg']/T['bp_sum']
        
        T = T[T['bp_sum']>0]
        
        S = T.groupby('pos').describe().unstack()
        col_names = [".".join(x) for x in S.columns]
        
        S.columns = col_names
        S['pos'] = S.index 
        
        if self.annot:
            S = self.add_annotation(S, T)
          
        sys.stderr.write("done\n")
        return S 

    def get_full_len_binned_stats(self, nbins=1000):
        sys.stderr.write("getting full len binned stats...")
        
        idxs, mu_by_bin = [], [] 
        for idx, g in self.grouped_by_gene:
            binned = stats.binned_statistic(np.arange(g['cvg'].shape[0]), g['cvg'], bins=nbins)
            idxs.append(np.ones(nbins)*idx)
            mu_by_bin.append(binned[0])
        
        sys.stderr.write("finished individual gene parsing...")

        n = self.full_stats.shape[0]
        T = pd.DataFrame({"idx":np.concatenate(idxs), 
                          "binned_mu":np.concatenate(mu_by_bin),
                          "pos":np.tile(np.arange(nbins)+1,n)})
        
        T = pd.merge(T, self.full_stats, left_on="idx", right_on="idx")
        T = T[T['mean']>0]
        
        mu_sum = pd.DataFrame(T.groupby("idx").apply(lambda x: np.sum(x['binned_mu'])))
        mu_sum['idx'] = mu_sum.index
        mu_sum.columns = ['mu_sum', "idx"]
        T = pd.merge(T,mu_sum, left_on="idx", right_on="idx")
        T['normalized'] = T['binned_mu']/T['mu_sum']
        
        T = T[T['mu_sum']>0]
         
        S = T.groupby('pos').describe().unstack()
        col_names = [".".join(x) for x in S.columns]
        
        S.columns = col_names
        S['pos'] = S.index 
        
        if self.annot:
            S = self.add_annotation(S, T)

        sys.stderr.write("done\n")
        return S 
 
if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_h5", required=True)
    parser.add_argument("--fn_out_full_len_meta", required=True)
    parser.add_argument("--fn_out_3p_bp_meta", required=True)
    parser.add_argument("--fn_out_stop_bp_meta", required=True)
    parser.add_argument("--fn_annotation", default=None)
    parser.add_argument("--fn_annotation_key", default=None)
    parser.add_argument("--name", required=True)
    o = parser.parse_args()
    
    h5_obj = h5FullGeneCvg(o.fn_h5, o.fn_annotation, o.fn_annotation_key)
    
    STOP_stats = h5_obj.get_stop_bp_meta()
    STOP_stats['name'] = o.name
    STOP_stats.to_csv(o.fn_out_stop_bp_meta, sep="\t", index=False, compression="gzip")
    
    UTR_3p_stats = h5_obj.get_3p_bp_meta()
    UTR_3p_stats['name'] = o.name
    UTR_3p_stats.to_csv(o.fn_out_3p_bp_meta, sep="\t", index=False, compression="gzip")

    percentile_stats = h5_obj.get_full_len_binned_stats()
    percentile_stats['name'] = o.name
    percentile_stats.to_csv(o.fn_out_full_len_meta, sep="\t", index=False, compression="gzip")
    
