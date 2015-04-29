import sys
from collections import defaultdict
from fuzzywuzzy import process
import pandas as pd
import numpy as np


class ID_dict_cls(object):

    def __init__(self, **kwargs):
        self.geneID_to_geneName = kwargs.get("geneID_to_geneName")
        self.transcriptID_to_geneName= kwargs.get("transcriptID_to_geneName") 
        self.geneName_to_transcriptIDs  = kwargs.get("geneName_to_transcriptIDs")
        self.geneName_to_geneIDs = kwargs.get("geneName_to_geneIDs")
        self.geneID_to_loc = kwargs.get("geneID_to_loc")
        self.transcriptID_to_loc = kwargs.get("transcriptID_to_loc")
        
def get_ID_dicts(fn_gff):
    
    geneID_to_geneName = {}
    geneID_to_loc = {}

    transcriptID_to_geneName = {}
    transcriptID_to_loc = {}

    geneName_to_transcriptIDs = defaultdict(list)
    geneName_to_geneIDs = defaultdict(list)
     
    for l in open(fn_gff):
        if l[0]!="#":
            rest,info = l.rstrip().rsplit(None,1)
            contig, source, type, rest = l.rstrip().split(None,3)

            if type == "transcript":
                start, end, o1, strand, o2, info = rest.split()
                start, end = int(start), int(end)
                idict = {d.split("=")[0]:d.split("=")[1] for d in info.split(";")}
                if "gene_name" in idict:
                    gene_name = idict['gene_name']
                else:
                    gene_name = idict['ID']
                
                geneID = idict['geneID']
                geneID_to_geneName[geneID] = gene_name
                transcriptID_to_geneName[idict['ID']] = gene_name 

                transcriptID_to_loc[idict['ID']] = tuple([contig,start,end])
                if not geneID in geneID_to_loc: 
                    geneID_to_loc[geneID] = tuple([contig, start, end])
                else:
                    old_tup = geneID_to_loc[geneID] 
                    geneID_to_loc[geneID] = tuple([contig, min(start,old_tup[1]), max(end,old_tup[2])])

                if not idict['ID'] in geneName_to_transcriptIDs[gene_name]:
                    geneName_to_transcriptIDs[gene_name].append(idict['ID'])

                if not geneID in geneName_to_geneIDs[gene_name]:
                    geneName_to_geneIDs[gene_name].append(geneID)
    
    return ID_dict_cls(geneID_to_geneName=geneID_to_geneName, 
                       transcriptID_to_geneName=transcriptID_to_geneName, 
                       geneName_to_transcriptIDs=geneName_to_transcriptIDs, 
                       geneName_to_geneIDs=geneName_to_geneIDs,
                       geneID_to_loc=geneID_to_loc,
                       transcriptID_to_loc=transcriptID_to_loc)

class RsemParser(object):

    def __init__(self,fn_gff, gene_fns, sample_names, sample_groups, alt_names=None):
        
        self.samples = sample_names
        self.sample_groups = sample_groups
        
        self.alt_names = None
        
        if not alt_names:
            self.alt_names = {s:s for s in self.samples}
        else:
            self.alt_names = {self.samples[i]:alt_names[i] for i in range(len(alt_names))}
        
        
        ID_dicts = get_ID_dicts(fn_gff)
        self.ID_dicts = ID_dicts

        self.geneID_to_geneName = ID_dicts.geneID_to_geneName
        self.transcriptID_to_geneName = ID_dicts.transcriptID_to_geneName
        self.geneName_to_transcriptIDs =  ID_dicts.geneName_to_transcriptIDs
        self.geneName_to_geneIDs =  ID_dicts.geneName_to_geneIDs
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
