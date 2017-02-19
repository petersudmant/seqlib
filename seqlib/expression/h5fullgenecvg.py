import tables
import numpy as np
import pandas as pd
import argparse
import sys
import pdb
import scipy.stats as stats
import logging
import pysam
from coveragedata import *
from gtf_to_genes import *

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
        
        self.a_TID = self.h5.createEArray(self.h5.root, 
                                          'transcriptID', 
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
        self.TIDs = []
        self.lengths = []
        self.starts = []
        self.stops = []
    
    def extend(self, cvg_obj):

        self.a_cvg.append(cvg_obj.RNA_cvg_view)
        self.a_idx.append(np.ones(cvg_obj.RNA_cvg_view.shape[0])*self.curr_idx)
        
        #self.GeneIDs.append(cvg_obj.g.gene_id)
        self.GeneIDs.append(cvg_obj.gene_id)
        self.TIDs.append(cvg_obj.TID)
        self.lengths.append(cvg_obj.RNA_cvg_view.shape[0])
        self.starts.append(cvg_obj.RNA_cvg_view_START)
        self.stops.append(cvg_obj.RNA_cvg_view_STOP)

        self.curr_idx+=1
    
    def close(self):
        
        self.a_GeneID.append(np.array(self.GeneIDs))
        self.a_TID.append(np.array(self.TIDs))
        self.a_length.append(np.array(self.lengths))
        self.a_start.append(np.array(self.starts))
        self.a_stop.append(np.array(self.stops))

        self.h5.close()

class h5FullGeneCvg(object):
    
    def __init__(self, fn, **kwargs):
        """
        idx = 1 idx per transcript / gene
        """
        
        fn_annots = kwargs.get("fn_annots")
        annot_key_vals = kwargs.get("annot_key_vals")

        sys.stderr.write("loading {fn}...".format(fn=fn))
        self.h5 = tables.openFile(fn, mode='r')
        
        self.cvg = self.h5.root.cvg[:]
        self.idx = self.h5.root.idx[:]

        self.GeneID = self.h5.root.geneID[:]
        self.transcriptID = self.h5.root.transcriptID[:]
        self.length = self.h5.root.length[:]
        self.start = self.h5.root.start[:]
        self.stop = self.h5.root.stop[:]
        
        self.inf = pd.DataFrame({"GeneID":self.GeneID,
                                 "transcriptID":self.transcriptID,
                                 "length":self.length,
                                 "start":self.start,
                                 "stop":self.stop, 
                                 "idx":np.unique(self.idx)})
        
        self.inf = self.inf.set_index("idx", drop=False)
        self.trancriptID_to_idx = dict(zip(self.inf.transcriptID,self.inf.index))
        
        """
        if "filtering" then filter everything here...
        """
   
        self.table = pd.DataFrame({"cvg":self.cvg, "idx":self.idx})
        
        self.grouped_by_gene = self.table.groupby('idx') 
        self.full_stats = self.grouped_by_gene.agg([np.sum, 
                                                    np.mean, 
                                                    np.std, 
                                                    np.median]).reset_index(level=2, drop=True)
        
        self.full_stats.columns = ["sum", "mean", "std", "median"]
        self.full_stats['idx'] = self.full_stats.index
        
        self.annot = False
        
        if fn_annots:
            self.annot = True
            for i, fn in enumerate(fn_annots):
                key,value = annot_key_vals[2*i], annot_key_vals[(2*i)+1] 
                t = pd.read_csv(fn, header=0, sep="\t")
                t_annot = pd.merge(self.inf, t, left_on = "GeneID", right_on=key, how="left")[['idx',value]]
                t_annot.columns = ['idx', "annot_%s"%value]
                self.full_stats = pd.merge(self.full_stats, t_annot)

        sys.stderr.write("done\n")



    def get_stop_bp_meta(self, n_bp=300, R_tc_filter = None):
        sys.stderr.write("getting stop codon centered...")
        
        idxs, cvg_by_bp, R_tcs = [], [], []
        n=0
        for idx, g in self.grouped_by_gene:
            
            stop = self.inf.ix[idx]['stop'] 
            length = self.inf.ix[idx]['length']
        
            mu_post_stop = np.mean(g.cvg[stop:])
            mu_pre_stop = np.mean(g.cvg[:stop])
            r_R_tc = max(1,mu_post_stop)/max(1,mu_pre_stop)
            R_tc  = np.log(max(1,mu_post_stop)/max(1,mu_pre_stop))/np.log(2)
            
            if stop >=n_bp and length >= stop+n_bp+1:
                cvg_by_bp.append(g.cvg[(stop-n_bp):(stop+n_bp+1)])
                idxs.append(np.ones(2*n_bp+1)*idx)
                R_tcs.append(np.ones(2*n_bp+1)*R_tc)
                n+=1
            
        sys.stderr.write("finished individual gene parsing...")
        
        T = pd.DataFrame({"idx":np.concatenate(idxs), 
                          "R_tc":np.concatenate(R_tcs),
                          "bp_cvg":np.concatenate(cvg_by_bp),
                          "pos":np.tile(np.arange(-n_bp,n_bp+1),n)})
        
        T = pd.merge(T, self.full_stats, left_on="idx", right_on="idx")
        T = T[T['mean']>1]

        if R_tc_filter is not None:
            T = T[T['R_tc'].apply(R_tc_filter)]
        
        bp_sum = pd.DataFrame(T.groupby("idx").apply(lambda x: np.sum(x['bp_cvg'])))
        bp_sum['idx'] = bp_sum.index
        bp_sum.columns = ['bp_sum', "idx"]
        T = pd.merge(T, bp_sum, left_on="idx", right_on="idx")
        T['bp_normalized'] = T['bp_cvg']/T['sum']
        
        T = T[T['bp_sum']>0]
        
        S = T.groupby('pos').describe().unstack()
        col_names = [".".join(x) for x in S.columns]
        
        S.columns = col_names
        S['pos'] = S.index 
        
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
        
          
        sys.stderr.write("done\n")
        return S 

    def get_CDS_UTR_split_binned_meta(self, CDS_bins=100, UTR_bins=100):
        sys.stderr.write("getting full len binned stats...")
        n_bins_total = CDS_bins + UTR_bins 
        #########
        #########
        #########
        #########

        idxs = [] 
        CDS_UTR_by_bin = []
        R_tcs = []
        
        for idx, g in self.grouped_by_gene:
            stop = self.inf.ix[idx]['stop'] 
            length = self.inf.ix[idx]['length']
            
            CDS_binned = stats.binned_statistic(np.arange(stop), g['cvg'][:stop], bins=CDS_bins)
            UTR_binned = stats.binned_statistic(np.arange(length-stop), g['cvg'][stop:], bins=UTR_bins)
            
            idxs.append(np.ones(n_bins_total)*idx)
            CDS_UTR_by_bin.append(CDS_binned[0])
            CDS_UTR_by_bin.append(UTR_binned[0])
            
            ###########
            
        sys.stderr.write("finished individual gene parsing...")

        n = self.full_stats.shape[0]
        pos=np.tile(np.arange(n_bins_total),n)
        T = pd.DataFrame({"idx":np.concatenate(idxs), 
                          "binned_cvg":np.concatenate(CDS_UTR_by_bin),
                          "pos":pos,
                          "pos_offset":pos-CDS_bins})
        
        T = pd.merge(T, self.full_stats, left_on="idx", right_on="idx")
        #T = T[T['mean']>0]
        
        cvg_sum = pd.DataFrame(T.groupby("idx").apply(lambda x: np.sum(x['binned_cvg'])))
        cvg_sum['idx'] = cvg_sum.index
        cvg_sum.columns = ['cvg_sum', "idx"]
        T = pd.merge(T,cvg_sum, left_on="idx", right_on="idx")
        T['normalized'] = T['binned_cvg']/T['cvg_sum']
        
        #T = T[T['cvg_sum']>0]
         
        S = T.groupby('pos').describe().unstack()
        col_names = [".".join(x) for x in S.columns]
        
        S.columns = col_names
        S['pos'] = S.index 
        

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
        

        sys.stderr.write("done\n")
        return S 
    
    def get_stop_bp_heatmap(self, nt=300):
        sys.stderr.write("getting stop codon centered...")
        
        idxs, cvg_by_bp, R_tcs = [], [], []
        n=0
        for idx, g in self.grouped_by_gene:
            if idx % 100==0:
                sys.stderr.write("%d..."%idx)

            stop = self.inf.ix[idx]['stop'] 
            length = self.inf.ix[idx]['length']
        
            mu_post_stop = np.mean(g.cvg[stop:])
            mu_pre_stop = np.mean(g.cvg[:stop])
            r_R_tc = max(1,mu_post_stop)/max(1,mu_pre_stop)
            R_tc  = np.log(max(1,mu_post_stop)/max(1,mu_pre_stop))/np.log(2)
            
            if stop >=nt and length >= stop+nt+1:
                cvg_by_bp.append(g.cvg[(stop-nt):(stop+nt+1)])
                idxs.append(np.ones(2*nt+1)*idx)
                R_tcs.append(np.ones(2*nt+1)*R_tc)
                n+=1
            
        sys.stderr.write("\nfinished individual gene parsing...")
        
        T = pd.DataFrame({"idx":np.concatenate(idxs), 
                          "R_tc":np.concatenate(R_tcs),
                          "bp_cvg":np.concatenate(cvg_by_bp),
                          "pos":np.tile(np.arange(-nt,nt+1),n)})
        
        T = pd.merge(T, self.full_stats, left_on="idx", right_on="idx")

        sys.stderr.write("done\n")
        return T
    
    def get_stop_binned_heatmap(self, bins=100):
        sys.stderr.write("getting STOP centered binned stats...")
        n_bins_total = 2*bins
        
        CDS_bins = bins
        UTR_bins = bins
        
        idxs = [] 
        CDS_UTR_by_bin = []
        R_tcs = []
        
        for idx, g in self.grouped_by_gene:
            if idx % 100==0:
                sys.stderr.write("%d..."%idx)
            
            stop = self.inf.ix[idx]['stop'] 
            length = self.inf.ix[idx]['length']


            mu_post_stop = np.mean(g.cvg[stop:])
            mu_pre_stop = np.mean(g.cvg[:stop])
            R_tc  = np.log(max(1,mu_post_stop)/max(1,mu_pre_stop))/np.log(2)
            
            R_tcs.append(np.ones(2*bins)*R_tc)
            CDS_binned = stats.binned_statistic(np.arange(stop), g['cvg'][:stop], bins=CDS_bins)
            UTR_binned = stats.binned_statistic(np.arange(length-stop), g['cvg'][stop:], bins=UTR_bins)
            
            idxs.append(np.ones(n_bins_total)*idx)
            CDS_UTR_by_bin.append(CDS_binned[0])
            CDS_UTR_by_bin.append(UTR_binned[0])
            
        sys.stderr.write("\nfinished individual gene parsing...")

        n = self.full_stats.shape[0]
        pos=np.tile(np.arange(n_bins_total),n)
        T = pd.DataFrame({"idx":np.concatenate(idxs), 
                          "binned_cvg":np.concatenate(CDS_UTR_by_bin),
                          "R_tc":np.concatenate(R_tcs),
                          "pos":pos,
                          "pos_offset":pos-CDS_bins})
        
        T = pd.merge(T, self.full_stats, left_on="idx", right_on="idx")

        sys.stderr.write("done\n")

        return T 


    def get_R_tc(self):

        sys.stderr.write("getting R centered on termination (STOP) codon centered...")
        
        outrows = []
        n=0
        for idx, g in self.grouped_by_gene:
            stop = self.inf.ix[idx]['stop'] 
            length = self.inf.ix[idx]['length']
            gene = self.inf.ix[idx]['GeneID']

            mu_post_stop = np.mean(g.cvg[stop:])
            mu_pre_stop = np.mean(g.cvg[:stop])
            med_post_stop = np.median(g.cvg[stop:])
            med_pre_stop= np.median(g.cvg[:stop])
            outrows.append({"idx":idx,
                            "gene":gene,
                            "length":length,
                            "mu_post_stop":mu_post_stop,
                            "mu_pre_stop":mu_pre_stop,
                            "med_post_stop":med_post_stop,
                            "med_pre_stop":med_pre_stop,
                            "R_tc":max(1,mu_post_stop)/max(1,mu_pre_stop),
                            "R_tc_median":max(1,med_post_stop)/max(1,med_pre_stop)})

        
        sys.stderr.write("finished individual gene parsing...")
        
        T = pd.DataFrame(outrows)
        T = pd.merge(T, self.full_stats, left_on="idx", right_on="idx")
        
        sys.stderr.write("done\n")
        return T 


def build_h5(args):
    
    logger = logging.getLogger(args.fn_logfile)
    bamfile = pysam.AlignmentFile(args.fn_bam, 'rb')
    
    contig_lengths = {x[0]:x[1] for x in zip(bamfile.references, bamfile.lengths)}
    
    species_id, gtf_path, genes = get_indexed_genes_for_identifier(args.fn_gtf_index,
                                                                   logger, 
                                                                   args.gtf_ID)
    if args.meta_type=="constitutive_single_stop":
        cvg_objs_by_contig = get_cvg_objs_by_contig(genes,"constitutive_single_stop")
    elif args.meta_type=="all_transcripts": 
        cvg_objs_by_contig = get_cvg_objs_by_contig(genes,"transcript")

    full_h5 = h5FullGeneCvg_writer(args.fn_out)

    max_t=0
    total_assessed = 0
    for contig, cvg_objs in cvg_objs_by_contig.items():
        if contig not in contig_lengths:
            continue
        print("current contig: %s\tassessed: %d"%(contig, total_assessed))
        t=time.time()
        cvg_recarray = pysamstats.load_nondel_coverage(bamfile, 
                                                      chrom=contig, 
                                                      start=0, 
                                                      end=contig_lengths[contig])
        print("time to load contig: %fs"%(time.time()-t))
        sys.stdout.flush()

        for i, cvg_obj in enumerate(cvg_objs):
            total_assessed +=1
            
            t = time.time()
            cvg_obj.get_cvg(cvg_recarray, bamfile)
            if time.time()-t>max_t: 
                max_t = time.time()-t
                print("t=%f, max_t=%f"%(time.time()-t,max_t))
                sys.stdout.flush()
            
            full_h5.extend(cvg_obj)
        
    full_h5.close()

def heatmap(args):
    
    h5_obj = h5FullGeneCvg(args.fn_h5, 
                           fn_annots = args.fn_annotations, 
                           annot_key_vals = args.annotation_key_values)
    
    if args.type =="stop_nt": 
        STOP_heatmap = h5_obj.get_stop_bp_heatmap(nt=args.size)
        STOP_heatmap['sample_name'] = args.sample_name
        STOP_heatmap.to_csv(args.fn_out, 
                            sep="\t", 
                            index=False, 
                            compression="gzip")
    if args.type =="stop_binned": 
        STOP_heatmap = h5_obj.get_stop_binned_heatmap(bins=args.size)
        STOP_heatmap['sample_name'] = args.sample_name
        STOP_heatmap.to_csv(args.fn_out, 
                            sep="\t", 
                            index=False, 
                            compression="gzip")
    else:
        assert True, "no such type %s"%args.type

def metaplot(args):
    
    h5_obj = h5FullGeneCvg(args.fn_h5, 
                           fn_annots = args.fn_annotations, 
                           annot_key_vals = args.annotation_key_values)
    
    if args.type =="full_len":
        percentile_stats = h5_obj.get_full_len_binned_stats()
        percentile_stats['sample_name'] = args.sample_name
        percentile_stats.to_csv(args.fn_out, sep="\t", index=False, compression="gzip")
    
    elif args.type =="3p_nt": 
        UTR_3p_stats = h5_obj.get_3p_bp_meta()
        UTR_3p_stats['sample_name'] = args.sample_name
        UTR_3p_stats.to_csv(args.fn_out, sep="\t", index=False, compression="gzip")
    
    elif args.type =="stop_nt": 
        if args.lambda_fun is not None:
            STOP_stats = h5_obj.get_stop_bp_meta(R_tc_filter = eval(args.lambda_fun))
        else:
            STOP_stats = h5_obj.get_stop_bp_meta()
        STOP_stats['sample_name'] = args.sample_name
        STOP_stats.to_csv(args.fn_out, sep="\t", index=False, compression="gzip")
    
    elif args.type =="binned": 
        stats = h5_obj.get_CDS_UTR_split_binned_meta()
        stats['sample_name'] = args.sample_name
        stats.to_csv(args.fn_out, sep="\t", index=False, compression="gzip")
    
def Rtc(args):
    
    h5_obj = h5FullGeneCvg(args.fn_h5, 
                           fn_annots = args.fn_annotations, 
                           annot_key_vals = args.annotation_key_values)

    R_tc = h5_obj.get_R_tc()
    R_tc['sample_name'] = args.sample_name
    R_tc.to_csv(args.fn_out, sep="\t", index=False, compression="gzip")

def plotCoverage(args):
    tables = []
    for i, sample in enumerate(args.samples):
        fn_h5 = args.fn_h5s[i]
        h5_obj = h5FullGeneCvg(fn_h5, 
                               fn_annots = args.fn_annotations, 
                               annot_key_vals = args.annotation_key_values)
        idx = h5_obj.trancriptID_to_idx[args.transcript]
        cvg = h5_obj.cvg[h5_obj.idx==idx]
        tid_inf = h5_obj.inf.ix[idx]
        
        T_out = pd.DataFrame({"cvg":cvg,
                              "pos":np.arange(cvg.shape[0])})
        T_out['GeneID'] = tid_inf['GeneID']
        T_out['transcriptID'] = tid_inf['transcriptID']
        T_out['length'] = tid_inf['length']
        T_out['start'] = tid_inf['start']
        T_out['stop'] = tid_inf['stop']

        T_out['sample'] = sample
        T_out['gene'] = args.gene
        tables.append(T_out)
    T = pd.concat(tables)
    T.to_csv(args.fn_out, sep="\t", index=False, compression="gzip")

if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers() 
    
    #create h5
    parser_create = subparsers.add_parser("build_h5")
    parser_create.add_argument("--fn_bam", required=True)
    parser_create.add_argument("--fn_out", required=True)
    parser_create.add_argument("--fn_gtf_index", required=True)
    parser_create.add_argument("--gtf_ID", required=True)
    parser_create.add_argument("--meta_type", required=True, 
                                              choices=["constitutive_single_stop",
                                                       "all_transcripts"])
    parser_create.add_argument("--fn_logfile", default="/dev/null")
    parser_create.set_defaults(func=build_h5)
    
    #output meta plot measurements
    parser_metaplot = subparsers.add_parser("metaplot")
    parser_metaplot.add_argument("--fn_h5", required=True)
    parser_metaplot.add_argument("--type", required=True, choices=["full_len",
                                                                   "3p_nt",
                                                                   "stop_nt",
                                                                   "binned"])
    parser_metaplot.add_argument("--fn_out", required=True)
    parser_metaplot.add_argument("--sample_name", required=True)
    parser_metaplot.add_argument("--lambda_fun", required=False)
    parser_metaplot.add_argument("--fn_annotations", default=None, nargs = "*")
    parser_metaplot.add_argument("--annotation_key_values", default=None, nargs="*")
    parser_metaplot.set_defaults(func=metaplot)
    
    #output heatmap data
    parser_metaplot = subparsers.add_parser("heatmap")
    parser_metaplot.add_argument("--fn_h5", required=True)
    parser_metaplot.add_argument("--type", required=True, choices=["stop_nt", "stop_binned"])
    parser_metaplot.add_argument("--size", required=False, type=int)
    parser_metaplot.add_argument("--fn_out", required=True)
    parser_metaplot.add_argument("--sample_name", required=True)
    parser_metaplot.add_argument("--lambda_fun", required=False)
    parser_metaplot.add_argument("--fn_annotations", default=None, nargs = "*")
    parser_metaplot.add_argument("--annotation_key_values", default=None, nargs="*")
    parser_metaplot.set_defaults(func=heatmap)

    #output Rtc metrics
    parser_Rtc = subparsers.add_parser("Rtc")
    parser_Rtc.add_argument("--fn_h5", required=True)
    parser_Rtc.add_argument("--fn_out", required=True)
    parser_Rtc.add_argument("--lambda_fun", required=False)
    parser_Rtc.add_argument("--sample_name", required=True)
    parser_Rtc.add_argument("--fn_annotations", default=None, nargs = "*")
    parser_Rtc.add_argument("--annotation_key_values", default=None, nargs="*")
    parser_Rtc.set_defaults(func=Rtc)
    
    #output gene coverage 
    parser_Rtc = subparsers.add_parser("plotCoverage")
    parser_Rtc.add_argument("--fn_h5s", required=True, nargs="+")
    parser_Rtc.add_argument("--fn_out", required=True)
    parser_Rtc.add_argument("--lambda_fun", required=False)
    parser_Rtc.add_argument("--samples", required=True, nargs="+")
    parser_Rtc.add_argument("--gene", required=True)
    parser_Rtc.add_argument("--transcript", required=True)
    parser_Rtc.add_argument("--fn_annotations", default=None, nargs = "*")
    parser_Rtc.add_argument("--annotation_key_values", default=None, nargs="*")
    parser_Rtc.set_defaults(func=plotCoverage)
    
    
    args = parser.parse_args()
    args.func(args)
     
    
    
