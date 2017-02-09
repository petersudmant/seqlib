"""
DESIGNED TO BE USED WITH /home/psudmant/code/RNAplfold_marvin
see /home/psudmant/projects/trap_development/analysis/3p_UTR_frags/assess_frags/structure 
for details
"""

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
        
        self.a_pos = self.h5.createEArray(self.h5.root, 
                                          'cvg', 
                                          tables.Float32Atom(), 
                                          shape=(0,), 
                                          filters=filt,
                                          expectedrows=65000000)
        
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
        
        if fn_annots:
            self.annot = True
            for i, fn in enumerate(fn_annots):
                key,value = annot_key_vals[2*i], annot_key_vals[(2*i)+1] 
                t = pd.read_csv(fn, header=0, sep="\t")
                t_annot = pd.merge(self.inf, t, left_on = "GeneID", right_on=key, how="left")[['idx',value]]
                t_annot.columns = ['idx', "annot_%s"%value]
                self.full_stats = pd.merge(self.full_stats, t_annot)

        sys.stderr.write("done\n")
    

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


if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers() 
    
    #create h5
    parser_create = subparsers.add_parser("build_h5")
    parser_create.add_argument("--fn_plfold", required=True)
    parser_create.add_argument("--fn_out", required=True)
    parser_create.add_argument("--fn_gtf_index", required=True)
    parser_create.add_argument("--gtf_ID", required=True)
    parser_create.add_argument("--fn_logfile", default="/dev/null")
    parser_create.set_defaults(func=build_h5)
    
    args = parser.parse_args()
    args.func(args)
     
    
    
