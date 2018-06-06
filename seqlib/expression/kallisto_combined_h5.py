import tables
import numpy as np
import pandas as pd
import argparse
import sys
import pdb
import scipy.stats as stats
import logging
import pysam
import json

class h5KallistoCvg_writer(object):
    """
    Kallisto CVG writer
    stores kallisto TPMs for multiple samples
    includes:
        list of all sample names
            *max sample name is 128 chars
        for each transcript:
            kallisto_id 
            transcript_id
            eff_length
        for each sample:
            TPM
                *TPM = est_counts/(s*hd_eff_length) 
                *s=np.sum(x.est_counts/x.eff_lengths)


    """
    
    def __init__(self, **kwargs):
        
        t_kallisto_annot = pd.read_csv(kwargs['fn_kallisto_idx_annot'], 
                                       names=["kallisto_id", "transcript_id"],
                                       header=None,
                                       sep=" ")
        
        n_transcripts = t_kallisto_annot.shape[0]

        eg_kallisto_quant = tables.openFile(kwargs['kallisto_abundance_fns'][0], mode='r')
        
        self.h5 = tables.openFile(args.fn_out, mode='w')
        filt = tables.Filters(complevel=5, complib='blosc')
        
        #transcript data
        self.a_transcript_id = self.h5.createCArray(self.h5.root, 
                                             'transcript_id', 
                                             tables.StringAtom(32),
                                             shape=(n_transcripts,), 
                                             filters=filt)

        self.a_kallisto_id = self.h5.createCArray(self.h5.root, 
                                             'kallisto_id', 
                                             tables.Int64Atom(),
                                             shape=(n_transcripts,), 
                                             filters=filt)
        
        self.a_eff_length = self.h5.createCArray(self.h5.root, 
                                             'eff_length', 
                                             tables.Float64Atom(),
                                             shape=(n_transcripts,), 
                                             filters=filt)
        
        self.a_length = self.h5.createCArray(self.h5.root, 
                                             'length', 
                                             tables.Int32Atom(),
                                             shape=(n_transcripts,), 
                                             filters=filt)
        

        self.a_transcript_id[:] = t_kallisto_annot["transcript_id"].values
        self.a_kallisto_id[:] = t_kallisto_annot["kallisto_id"].values
        self.a_eff_length[:] = eg_kallisto_quant.root.aux.eff_lengths[:]
        self.a_length[:] = eg_kallisto_quant.root.aux.lengths[:]

        #samples - extendable
        self.a_sample = self.h5.createEArray(self.h5.root, 
                                             'sample', 
                                             tables.StringAtom(128),
                                             shape=(0,), 
                                             filters=filt,
                                             expectedrows=10000)


        self.a_sample.append(kwargs['samples'])
        
        #per sample quantifications - extendable


        self.a_TPM = self.h5.createEArray(self.h5.root, 
                                          'TPM', 
                                          tables.Float64Atom(), 
                                          shape=(0,n_transcripts), 
                                          filters=filt)
        
        
        
        for i, sample in enumerate(kwargs['samples']):
            if i%100==0:
                sys.stdout.write("{i}/{t}...".format(i=i,t=len(kwargs['samples'])))
                sys.stdout.flush()
            h5_kquant = tables.openFile(kwargs['kallisto_abundance_fns'][i], mode='r')
            s = np.sum(h5_kquant.root.est_counts[:]/h5_kquant.root.aux.eff_lengths[:])
            TPM = h5_kquant.root.est_counts[:]/(s*h5_kquant.root.aux.eff_lengths[:])*1e6
            self.h5.root.TPM.append(TPM.reshape((1,-1)))

    def extend(self, cvg_obj):
        """
            not implemented
            adds to already built file
        """
        return
    
    def close(self):
        self.h5.close()


def build_h5(args):
    #parser_create.add_argument("--fn_kallisto_idx_annot", required=True)

    samples = json.load(open(args.sample_json))['samples']
    kallisto_abundance_fns = ["{path}/{s}/abundance.h5".format(path=args.kallisto_abundance_path, s=s) for s in samples]
    h5 = h5KallistoCvg_writer(fn_kallisto_idx_annot = args.fn_kallisto_idx_annot, 
                              kallisto_abundance_fns = kallisto_abundance_fns,
                              fn_out=args.fn_out, 
                              samples=samples)
    h5.close()
   

def summarize_introns(args):
    
    h5 = tables.openFile(args.fn_h5, mode='r')
      

if __name__=="__main__":
    
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers() 
    
    #create h5
    parser_create = subparsers.add_parser("build_h5")
    parser_create.add_argument("--fn_out", required=True)
    parser_create.add_argument("--fn_kallisto_idx_annot", required=True)
    parser_create.add_argument("--kallisto_abundance_path", required=True)
    parser_create.add_argument("--sample_json", required=True)
    parser_create.set_defaults(func=build_h5)
    
    #create h5
    #parser_summarize_introns = subparsers.add_parser("summarize_introns")
    #parser_create.add_argument("--fn_out", required=True)
    #parser_create.add_argument("--fn_h5", required=True)
    #parser_create.set_defaults(func=summarize_introns)
    
    args = parser.parse_args()
    args.func(args)
