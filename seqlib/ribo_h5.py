import tables
import logging 
import numpy as np
import pandas as pd
import argparse
import sys
import pdb
import pysam
import scipy.stats as stats

from expression.coveragedata import *
from gtf_to_genes import *

class h5_ribo_writer(object):
    
    def __init__(self, fn):
        self.h5 = tables.openFile(fn, mode='w')
        
        filt = tables.Filters(complevel=5, complib='blosc')
        """
        each of these is one entry per read
        """
        self.pos = self.h5.createEArray(self.h5.root, 
                                        'pos', 
                                        tables.UInt32Atom(), 
                                        shape=(0,), 
                                        filters=filt,
                                        expectedrows=65000000)
        
        self.strand = self.h5.createEArray(self.h5.root, 
                                           'strand', 
                                           tables.Int8Atom(), 
                                           shape=(0,), 
                                           filters=filt,
                                           expectedrows=65000000)
        
        self.contig_idx = self.h5.createEArray(self.h5.root, 
                                               'contig_idx', 
                                               tables.UInt8Atom(), 
                                               shape=(0,), 
                                               filters=filt,
                                               expectedrows=65000000)
        
        self.length = self.h5.createEArray(self.h5.root, 
                                           'length', 
                                           tables.UInt8Atom(), 
                                           shape=(0,), 
                                           filters=filt,
                                           expectedrows=65000000)
        """
        meta data storage for contigs
        """

        self.contigs = self.h5.createEArray(self.h5.root, 
                                            'contigs', 
                                            tables.StringAtom(50),
                                            shape=(0,), 
                                            filters=filt,
                                            expectedrows=100)
        
        self.curr_idx = 0
        self.contig_list = []
    
    def populate(self, bamfile):
    
        contigs = {x[0]:x[1] for x in zip(bamfile.references, bamfile.lengths)}
        
        n = 0
        for contig, clen in contigs.items():
            if "chrUn" in contig: continue
            if "random" in contig: continue
                
            sys.stderr.write("loading {contig}...".format(contig=contig))
            strands = []
            lens = []
            poses = []

            for r in bamfile.fetch(contig):
                strand = not(r.is_reverse)
                pos = r.pos
                
                poses.append(pos)
                strands.append(strand)
                lens.append(r.query_alignment_length)

                n+=1
            
            l = len(poses)
            self.pos.append(np.array(poses,dtype='uint32'))
            self.strand.append(np.array(strands,dtype='int8'))
            self.length.append(np.array(lens,dtype='uint8'))
            self.contig_idx.append(np.ones(l)*self.curr_idx)
            
            self.contig_list.append(contig) 
            self.curr_idx+=1

    def close(self):
        self.contigs.append(np.array(self.contig_list))
        self.h5.close()

class h5_ribo(object):
    
    def __init__(self, fn, **kwargs):
        
        sys.stderr.write("loading {fn}...".format(fn=fn))
        self.h5 = tables.openFile(fn, mode='r')

        self.pos = self.h5.root.pos[:]
        self.strand = self.h5.root.strand[:]
        self.length = self.h5.root.length[:]
        
        self.contig_idx = self.h5.root.contig_idx[:]
        self.contigs = self.h5.root.contigs[:]
        
        sys.stderr.write("done\n")
    
    def get_RPF_counts(self, offsets):
        
        offset_poses = np.zeros(self.pos.shape[0], dtype='uint32')
        for size, offset in offsets.items(): 
            curr_size_idx_p = np.where((self.length==size)&(self.strand==1))
            curr_size_idx_m = np.where((self.length==size)&(self.strand==0))
            offset_poses[curr_size_idx_p] = self.pos[curr_size_idx_p]+offset+1
            offset_poses[curr_size_idx_m] = self.pos[curr_size_idx_m]+self.length[curr_size_idx_m]-offset

        non_zero_idx = np.where(offset_poses!=0)
        offset_poses = offset_poses[non_zero_idx]
        strand = self.strand[non_zero_idx]

        contig_idx = self.contig_idx[non_zero_idx]
        t = pd.DataFrame({"pos":offset_poses,
                          "strand":strand,
                          "contig_idx":contig_idx})

        t_grp = t.groupby(["pos","strand","contig_idx"])

        RPF_counts = t_grp.pos.size().to_frame().reset_index()
        RPF_counts.columns = ["pos", "strand", "contig_idx", "count"]
        RPF_counts = RPF_counts.sort_values(by=["contig_idx","pos"])
        RPF_counts['contig'] = RPF_counts.contig_idx.apply(lambda x: self.contigs[x])

        return RPF_counts

def create(args):
    bamfile = pysam.AlignmentFile(args.fn_bam, 'rb')
    new_h5 = h5_ribo_writer(args.fn_out)
    new_h5.populate(bamfile)
    new_h5.close()


def get_aSiteOffsets(aSiteOffset_list):
    
    aSiteOffsets = {} 
    for tup in aSiteOffset_list:
        size, offset = [int(x) for x in tup.split(":")]
        aSiteOffsets[size] = offset
    
    return aSiteOffsets


def get_cvg(contig, exons, strand, t_RPF):
    
    sum = 0
    l = 0
    for ex in exons:
        s,e = ex
        l+=e-s
        sum += np.sum(t_RPF[(t_RPF.contig==contig)&(t_RPF.pos>=s)&(t_RPF.pos<e)]['count'])

    return sum, l

def makeSummary(args):
    h5 = h5_ribo(args.fn_h5)
    offsets =  get_aSiteOffsets(args.aSiteOffsets)
    RPF_counts = h5.get_RPF_counts(offsets)
    
    logger = logging.getLogger(args.fn_logfile)
    s_id, path, genes = get_indexed_genes_for_identifier(args.fn_gtf_index,
                                                                   logger, 
                                                                   args.gtf_ID)
    cvg_objs_by_contig = get_cvg_objs_by_contig(genes,
                                                "constitutive_single_stop")
    outrows = []
    for contig, cvg_objs in cvg_objs_by_contig.items():
        for cvg_ob in cvg_objs:
            cding_cvg, cding_len = get_cvg(contig, cvg_ob.coding_exons, cvg_ob.strand, RPF_counts)
            
            outrows.append(
            pdb.set_trace()

def makeWig(args):
    
    h5 = h5_ribo(args.fn_h5)
    offsets =  get_aSiteOffsets(args.aSiteOffsets)
    color_by_strand = {0:"200,0,0",
                       1:"0,200,0"}
    
    RPF_counts = h5.get_RPF_counts(offsets)
    
    name = args.fn_out.split("/")[-1].split(".wig")[0]
    Fout = open(args.fn_out,'w')
    
    for strand in [0,1]:
        RPF_at_strand = RPF_counts[RPF_counts['strand']==strand]
        
        Fout.write("""track type=wiggle_0 name="{name}_strand{strand}" """
                   """visibility=full autoScale=on alwaysZero=on """
                   """color={color} priority=10\n""".format(name=name,
                                                            strand=strand,
                                                            color=color_by_strand[strand]))
        for i, contig in enumerate(h5.contigs):
            Fout.write("variableStep chrom={contig}\n".format(contig=contig))
            t_out = RPF_at_strand[RPF_at_strand.contig_idx==i][['pos','count']]
            t_out.to_csv(Fout, sep=" ", header=False, index=False)

    Fout.close()

if __name__=="__main__":

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers() 
    
    #create h5
    parser_create = subparsers.add_parser("create")
    parser_create.add_argument("--fn_bam", required=True)
    parser_create.add_argument("--fn_out", required=True)
    parser_create.set_defaults(func=create)
    
    #make wig from h5
    parser_makeWig = subparsers.add_parser("makeWig")
    parser_makeWig.add_argument("--fn_h5", required=True)
    parser_makeWig.add_argument("--fn_out", required=True)
    parser_makeWig.add_argument("--aSiteOffsets", required=True, nargs="+")
    parser_makeWig.set_defaults(func=makeWig)
    
    #make summary of coverage details
    parser_makeSum = subparsers.add_parser("makeSummary")
    parser_makeSum.add_argument("--fn_h5", required=True)
    parser_makeSum.add_argument("--fn_out", required=True)
    parser_makeSum.add_argument("--aSiteOffsets", required=True, nargs="+")
    parser_makeSum.add_argument("--fn_gtf_index", required=True)
    parser_makeSum.add_argument("--gtf_ID", required=True)
    parser_makeSum.add_argument("--fn_logfile", default='/dev/stderr')
    parser_makeSum.set_defaults(func=makeSummary)
    
    args = parser.parse_args()
    args.func(args)


