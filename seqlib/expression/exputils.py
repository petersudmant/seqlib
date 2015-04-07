import sys
from collections import defaultdict
from fuzzywuzzy import process
import pandas as pd
import numpy as np

class RsemParser(object):

    def __init__(self,fn_gff, gene_fns, sample_names, sample_groups, alt_names=None):
        
        self.samples = sample_names
        self.sample_groups = sample_groups
        
        self.alt_names = None
        
        if not alt_names:
            self.alt_names = {s:s for s in self.samples}
        else:
            self.alt_names = {self.samples[i]:alt_names[i] for i in range(len(alt_names))}
        
        geneID_to_geneName = {}
        transcriptID_to_geneName = {}
        geneName_to_transcriptIDs = defaultdict(list)
        geneName_to_geneIDs = defaultdict(list)
         
        for l in open(fn_gff):
            if l[0]!="#":
                rest,info = l.rstrip().rsplit(None,1)
                if ";" in info:
                    idict = {d.split("=")[0]:d.split("=")[1] for d in info.split(";")}
                    if "gene_name" in idict:
                        geneID_to_geneName[idict['geneID']] = idict['gene_name']
                        transcriptID_to_geneName[idict['ID']] = idict['gene_name']
                        if not idict['ID'] in geneName_to_transcriptIDs[idict['gene_name']]:
                            geneName_to_transcriptIDs[idict['gene_name']].append(idict['ID'])
                        if not idict['geneID'] in geneName_to_geneIDs[idict['gene_name']]:
                            geneName_to_geneIDs[idict['gene_name']].append(idict['geneID'])
        
        print("%d genes parsed"%(len(geneID_to_geneName.keys())), file=sys.stderr)

        self.geneID_to_geneName = geneID_to_geneName
        self.transcriptID_to_geneName = transcriptID_to_geneName
        self.geneName_to_transcriptIDs =  geneName_to_transcriptIDs
        self.geneName_to_geneIDs =  geneName_to_geneIDs
        self.geneNames = list(self.geneName_to_geneIDs.keys())

        gene_fpkms=None
        
        for i, fn in enumerate(gene_fns):
            t=pd.read_csv(fn,header=0,sep="\t")
            t=t[["gene_id", 'FPKM']]
            t.columns = ["gene_id", sample_names[i]]
            if i==0:
                gene_fpkms = t
            else:
                gene_fpkms = gene_fpkms.merge(t,on="gene_id")
        
        self.gene_fpkms = gene_fpkms.set_index(["gene_id"])

    
    def get_fuzzy_matching_genes(self, name, cutoff=90):
        
        fuzzy_matches = []
        for ftup in process.extract(name, self.geneNames):
            if ftup[1]>=cutoff:
                fuzzy_matches.append(ftup[0])
        return fuzzy_matches
    
    def get_gene_fpkms(self, genes, filter_to_max=True):
        """
        filter to max filters to one entry per gene name by the MAX mean
        across all samples - a hack, but, it works
        """
         
        fpkm_by_sample_by_gene = {}
        rows = []

        for g in genes:
            curr_mx = -1
            curr_fpkms = {}
            for ID in self.geneName_to_geneIDs[g]:
                fpkms = self.gene_fpkms.loc[ID]
                if filter_to_max: 
                    if np.mean(fpkms.values)>curr_mx:
                        curr_fpkms = {tuple([g,ID]):fpkms}
                        curr_mx = np.mean(fpkms.values)
                else: 
                    curr_fpkms[tuple([g,ID])] = fpkms

            for name_ID_tup, fpkm in curr_fpkms.items():
                for s in self.samples:
                    rows.append({"GeneName":name_ID_tup[0],
                                 "GeneID":name_ID_tup[1],
                                 "sample":s,
                                 "alt_name":self.alt_names[s],
                                 "group":self.sample_groups[s],
                                 "fpkm":fpkm[s]})

        return pd.DataFrame(rows)
