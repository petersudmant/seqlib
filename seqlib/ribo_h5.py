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
        """
        stores the position of the 5' end of each read
        """
        contigs = {x[0]:x[1] for x in zip(bamfile.references, bamfile.lengths)}
         
        n = 0
        for contig, clen in contigs.items():
            if "chrUn" in contig: continue
            if "random" in contig: continue
            
            sys.stderr.write("loading {contig}...".format(contig=contig))
            strands = []
            lens = []
            poses = []

            #for r in bamfile.fetch("chr5",137891642,137916204):
            for r in bamfile.fetch(contig):
                strand = not(r.is_reverse)
                pos = r.pos
                #NOW set pos to 5' end
                if strand:
                    poses.append(pos)
                else:
                    poses.append(pos+r.query_alignment_length)

                strands.append(strand)
                lens.append(r.query_alignment_length)
                #print(r)
                #print(pos, strand, r.query_alignment_length, r.cigartuples)
                #pdb.set_trace()
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
        self.offsets = None        
 
        self.contig_to_idx = {}
        for i in range(self.contigs.shape[0]):
            self.contig_to_idx[self.contigs[i]] = i

        sys.stderr.write("done\n")
         
    def load_offsets(self, fn):
        t = pd.read_csv(fn, sep=" ", header=0)
        self.offsets = {}
        for i, row in t.iterrows():
            assert row['manual'] != 99
            self.offsets[row['len']] = row['manual']
             
    def get_offset_contig_arrays(self, contig):
        assert self.offsets is not None, "need to assign offsets"        
        
        contig_pos = self.pos[self.contig_idx == self.contig_to_idx[contig]]
        contig_strand = self.strand[self.contig_idx == self.contig_to_idx[contig]]
        contig_length = self.length[self.contig_idx == self.contig_to_idx[contig]]
        
        offset_contig_pos = np.zeros(contig_pos.shape[0],dtype='uint32')
        for l, o in self.offsets.items():
            w_pos = np.where((contig_length==l)&(contig_strand==1))
            w_neg = np.where((contig_length==l)&(contig_strand==0))
            offset_contig_pos[w_pos] = contig_pos[w_pos]+o
            offset_contig_pos[w_neg] = contig_pos[w_neg]-o
        
        #ONLY KEEP THOSE THAT are amongst the lengths we want 
        #!!! need to assign strand before I mutate offset_contig_pos!
        contig_strand = contig_strand[offset_contig_pos!=0]
        offset_contig_pos = offset_contig_pos[offset_contig_pos!=0]
        
        s_neg = np.where(contig_strand==0)
        s_pos = np.where(contig_strand==1)
        
        return {0:offset_contig_pos[s_neg],
                1:offset_contig_pos[s_pos]}
    
    def get_offset_counts_by_contig(self, contig):
        pos_arrays = self.get_offset_contig_arrays(contig)
        
        count_arrays = {}
        for strand, pos_array in pos_arrays.items():        
            uniq, counts = np.unique(pos_array, return_counts = True)
            count_arrays[strand] = {"pos":uniq,
                                    "count":counts}
        return count_arrays
    

    def get_contig_arrays(self, contig):

        contig_pos = self.pos[self.contig_idx == self.contig_to_idx[contig]]
        contig_strand = self.strand[self.contig_idx == self.contig_to_idx[contig]]
        contig_length = self.length[self.contig_idx == self.contig_to_idx[contig]]

        s_neg = np.where(contig_strand==0)
        s_pos = np.where(contig_strand==1)
        
        return {0:
                {   "pos":contig_pos[s_neg],
                    "length":contig_length[s_neg]
                },
               1:
                {   "pos":contig_pos[s_pos],
                    "length":contig_length[s_pos]
                }
               }

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

def calibrate(args, width=50):
    h5 = h5_ribo(args.fn_h5)
    
    sys.stderr.write("loading gene annotations...")
    logger = logging.getLogger(args.fn_logfile)
    s_id, path, genes = get_indexed_genes_for_identifier(args.fn_gtf_index,
                                                                   logger, 
                                                                   args.gtf_ID)
    sys.stderr.write("done\n")
    cvg_objs_by_contig = get_cvg_objs_by_contig(genes,
                                                "transcript")
    outrows = []
    
    stop_c_by_size = {}
    start_c_by_size = {}
    #counts_by_type_len = {"start":{},"stop":{}} 
    counts_by_type_len = {"start":{0:{},1:{}},"stop":{0:{},1:{}}}

    for contig, cvg_objs in cvg_objs_by_contig.items():
        if not contig in h5.contig_to_idx: continue
         
        reads_by_strand = h5.get_contig_arrays(contig)

        for cvg_ob in cvg_objs:
            if cvg_ob.strand == 1:
                start, stop = cvg_ob.coding_exons[0][0], cvg_ob.coding_exons[-1][1] 
            else:
                start, stop = cvg_ob.coding_exons[-1][1], cvg_ob.coding_exons[0][0]
            
            all_pos = reads_by_strand[cvg_ob.strand]['pos']
            for typ,mid in {"start":start, "stop":stop}.items():
                w = np.where((all_pos>=mid-width)&(all_pos<=mid+width))
                poses = all_pos[w]
                lens = reads_by_strand[cvg_ob.strand]['length'][w]
                for l in np.unique(lens):
                    if not l in counts_by_type_len[typ][cvg_ob.strand]:
                        counts_by_type_len[typ][cvg_ob.strand][l]= np.zeros(2*width+1)
                    uniq, counts = np.unique(poses[lens==l], return_counts = True)
                    uniq = uniq.astype("int64") #so we can have -vs 
                    if cvg_ob.strand==1: 
                        centered_poses = uniq-mid
                    else:
                        centered_poses = mid-uniq
                    counts_by_type_len[typ][cvg_ob.strand][l][centered_poses+width]+= counts
    
    outrows = []
    for typ, counts_by_len_s in counts_by_type_len.items():
        for strand, counts_by_len in counts_by_len_s.items():
            for len, counts in counts_by_len.items():
                for i, pos in enumerate(range(-width,width)):
                    outrows.append({"type":typ,
                                    "len":len,
                                    "pos": pos,
                                    "strand":strand,
                                    "count":counts[i]})
                           
    t = pd.DataFrame(outrows)
    t.to_csv(args.fn_out, sep="\t", index=False)


def get_cvg(exons, strand, counts_by_contig):
    l = 0
    cvg = 0
    for e in exons:
        s,e = e
        l+=e-s
        w = np.where((counts_by_contig[strand]['pos']>=s)&(counts_by_contig[strand]['pos']<e))
        cvg += np.sum(counts_by_contig[strand]['count'][w])
    return cvg, l

def get_coverage_info(counts_by_contig, cvg_ob):
    
    if cvg_ob.strand == 1:
        start, stop = cvg_ob.coding_exons[0][0], cvg_ob.coding_exons[-1][1] 
        start_e, stop_e = [[start,start+3]],[[stop-2, stop+1]]
    else:
        start, stop = cvg_ob.coding_exons[-1][1], cvg_ob.coding_exons[0][0]
        start_e, stop_e = [[start-2,start+1]], [[stop, stop+3]]
    
    strand = cvg_ob.strand 
    CDS_cvg, CDS_l = get_cvg(cvg_ob.coding_exons, strand, counts_by_contig)
    UTR_3p_cvg, UTR_3p_l = get_cvg(cvg_ob.UTR_3p_exons, strand, counts_by_contig)
    UTR_5p_cvg, UTR_5p_l = get_cvg(cvg_ob.UTR_5p_exons, strand, counts_by_contig)
    STOP_cvg, l = get_cvg(stop_e, strand, counts_by_contig)
    START_cvg, l = get_cvg(start_e, strand, counts_by_contig)
    
    return {"CDS_cvg":CDS_cvg,
            "CDS_len":CDS_l,
            "UTR_5p_cvg":UTR_5p_cvg,
            "UTR_5p_len":UTR_5p_l,
            "UTR_3p_cvg":UTR_3p_cvg,
            "UTR_3p_len":UTR_3p_l,
            "START_cvg":START_cvg,
            "STOP_cvg":STOP_cvg}

def get_bp_coverage(counts_by_contig, cvg_ob, width, feature_pos):
    
    strand = cvg_ob.strand
    s,e = feature_pos-width, feature_pos+width+1
    cvg_vect = np.zeros(2*width+1)
    
    w = np.where((counts_by_contig[strand]['pos']>=s)&(counts_by_contig[strand]['pos']<e))
    poses = counts_by_contig[strand]['pos'][w]
    cvg = counts_by_contig[strand]['count'][w]

    if strand==1: 
        centered_poses = poses-feature_pos
    else:
        centered_poses = feature_pos-poses
    
    cvg_vect[centered_poses+width] = cvg
    
    return cvg_vect
    
    

def RPFCountTable(args):
    h5 = h5_ribo(args.fn_h5)
    h5.load_offsets(args.fn_aSiteOffsets)
    
    gene_id_subset = None 
    if args.fn_gene_subset is not None:
        t_subset = pd.read_csv(args.fn_gene_subset, header=0, sep="\t")
        gene_id_subset = t_subset.gene_id.values
 
    sys.stderr.write("loading gene annotations...")
    logger = logging.getLogger(args.fn_logfile)
    s_id, path, genes = get_indexed_genes_for_identifier(args.fn_gtf_index,
                                                                   logger, 
                                                                   args.gtf_ID)
    sys.stderr.write("done\n")
    cvg_objs_by_contig = get_cvg_objs_by_contig(genes,
                                                "transcript")
    outrows = []
    for contig, cvg_objs in cvg_objs_by_contig.items():
        if not contig in h5.contig_to_idx: continue
        
        sys.stderr.write("{contig}...".format(contig=contig))
        counts_by_contig = h5.get_offset_counts_by_contig(contig) #pos #count
        
        for cvg_ob in cvg_objs:
            if gene_id_subset is not None:
                if not cvg_ob.gene_id in gene_id_subset: continue

            pos = -1
            if args.feature=="STOP":
                if cvg_ob.strand == 1:
                    pos = cvg_ob.coding_exons[-1][1] 
                else:
                    pos = cvg_ob.coding_exons[0][0]
            elif args.feature=="START":
                if cvg_ob.strand == 1:
                    pos = cvg_ob.coding_exons[0][0]
                else:
                    pos = cvg_ob.coding_exons[-1][1]
            elif args.feature=="POLYA":
                if cvg_ob.strand == 1:
                    pos = cvg_ob.exons[-1][1]
                else:
                    pos = cvg_ob.exons[0][0]
            else:
                assert False, "no method for feature: %s"%(args.feature) 

            cvg_vect =  get_bp_coverage(counts_by_contig, cvg_ob, args.width, pos)
            cvg_info = get_coverage_info(counts_by_contig, cvg_ob)

            for i, pos in enumerate(range(-args.width,args.width)):
                outrows.append({"tid":cvg_ob.TID,
                                "gene_id":cvg_ob.gene_id,
                                "gene_name":cvg_ob.g.names[0],
                                "cvg":cvg_vect[i],
                                "CDS_cvg":cvg_info['CDS_cvg'],
                                "CDS_len":cvg_info['CDS_len'],
                                "pos":pos,
                                "feature":args.feature})

    t = pd.DataFrame(outrows)
    t.to_csv(args.fn_out, sep="\t", index=False)

def makeSummary(args):

    h5 = h5_ribo(args.fn_h5)
    h5.load_offsets(args.fn_aSiteOffsets)
    
    sys.stderr.write("loading gene annotations...")
    logger = logging.getLogger(args.fn_logfile)
    s_id, path, genes = get_indexed_genes_for_identifier(args.fn_gtf_index,
                                                                   logger, 
                                                                   args.gtf_ID)
    sys.stderr.write("done\n")
    cvg_objs_by_contig = get_cvg_objs_by_contig(genes,
                                                "transcript")

    outrows = []
    for contig, cvg_objs in cvg_objs_by_contig.items():
        if not contig in h5.contig_to_idx: continue
        
        sys.stderr.write("{contig}...".format(contig=contig))
        counts_by_contig = h5.get_offset_counts_by_contig(contig) #pos #count
        
        for cvg_ob in cvg_objs:
            cvg_inf = get_coverage_info(counts_by_contig, cvg_ob)
            cvg_inf.update({"tid":cvg_ob.TID,
                            "gene_id":cvg_ob.gene_id,
                            "gene_name":cvg_ob.g.names[0]})
            outrows.append(cvg_inf)
    t = pd.DataFrame(outrows)
    t.to_csv(args.fn_out, sep="\t", index=False)


def makeWig(args):
    
    h5 = h5_ribo(args.fn_h5)
    h5.load_offsets(args.fn_aSiteOffsets)
    
    color_by_strand = {0:"200,0,0",
                       1:"0,200,0"}
    
    name = args.fn_out.split("/")[-1].split(".wig")[0]
    Fout = open(args.fn_out,'w')
    
    
 
    outtables_by_strand_contig = {0:{},
                                  1:{}}
     
    for i, contig in enumerate(h5.contigs):
        counts_by_strand = h5.get_offset_counts_by_contig(contig)
        for strand, count_inf in counts_by_strand.items():
            t_out = pd.DataFrame({"pos":count_inf['pos'], "count":count_inf['count']})
            outtables_by_strand_contig[strand][contig] = t_out
    
    for strand, tables_by_contig in outtables_by_strand_contig.items():
        Fout.write("""track type=wiggle_0 name="{name}_strand{strand}" """
                   """visibility=full autoScale=on alwaysZero=on """
                   """color={color} priority=10\n""".format(name=name,
                                                            strand=strand,
                                                            color=color_by_strand[strand]))
        for contig, t in tables_by_contig.items():
            Fout.write("variableStep chrom={contig}\n".format(contig=contig))
            t.to_csv(Fout, sep=" ", header=False, index=False, columns=["pos","count"])
            
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
    parser_makeWig.add_argument("--fn_aSiteOffsets", required=True)
    parser_makeWig.set_defaults(func=makeWig)

    #callibrate
    parser_makeSum = subparsers.add_parser("calibrate")
    parser_makeSum.add_argument("--fn_h5", required=True)
    parser_makeSum.add_argument("--fn_out", required=True)
    parser_makeSum.add_argument("--fn_gtf_index", required=True)
    parser_makeSum.add_argument("--gtf_ID", required=True)
    parser_makeSum.add_argument("--fn_logfile", default='/dev/stderr')
    parser_makeSum.set_defaults(func=calibrate)

    #make summary of coverage details
    parser_makeSum = subparsers.add_parser("makeSummary")
    parser_makeSum.add_argument("--fn_h5", required=True)
    parser_makeSum.add_argument("--fn_out", required=True)
    parser_makeSum.add_argument("--fn_aSiteOffsets", required=True)
    parser_makeSum.add_argument("--fn_gtf_index", required=True)
    parser_makeSum.add_argument("--gtf_ID", required=True)
    parser_makeSum.add_argument("--fn_logfile", default='/dev/stderr')
    parser_makeSum.set_defaults(func=makeSummary)
    
    #output feature centered counts
    parser_makeSum = subparsers.add_parser("RPFCountTable")
    parser_makeSum.add_argument("--fn_h5", required=True)
    parser_makeSum.add_argument("--fn_out", required=True)
    parser_makeSum.add_argument("--fn_aSiteOffsets", required=True)
    parser_makeSum.add_argument("--fn_gtf_index", required=True)
    parser_makeSum.add_argument("--gtf_ID", required=True)
    parser_makeSum.add_argument("--fn_logfile", default='/dev/stderr')
    parser_makeSum.add_argument("--feature", required=True, choices = ["START","STOP","POLYA"])
    parser_makeSum.add_argument("--width", default=50, type=int)
    parser_makeSum.add_argument("--fn_gene_subset", required=False, default=None)
    parser_makeSum.set_defaults(func=RPFCountTable)
    
    args = parser.parse_args()
    args.func(args)


