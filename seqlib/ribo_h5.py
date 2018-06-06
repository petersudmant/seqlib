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
import timeit

class h5_ribo_writer(object):
    
    def __init__(self, fn):
        self.h5 = tables.openFile(fn, mode='w')
        self.fn = fn        
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
   
    def load_se(self, contig, clen, bamfile):
        strands = []
        lens = []
        poses = []
        
        for r in bamfile.fetch(contig):
            strand = not(r.is_reverse)
            pos = r.pos
            qual = r.mapping_quality
            #only take unique reads - STAR -> 255 for uniq
            if qual!=255: continue
            
            #if only using non-softclipped reads, check
            if nosoftclipping:
                if len(r.cigartuples)>1:
                    continue

            #NOW set pos to 5' end
            if strand:
                poses.append(pos)
            else:
                poses.append(pos+r.query_alignment_length)

            strands.append(strand)
            lens.append(r.query_alignment_length)
        
        return strands, lens, poses

    def load_pe(self, contig, clen, bamfile, pe_3p):
        """
        for pe reads, we take read1 5' end and read2 3' end
        note for pe reads the lengths are only useful to get the proper end
        because they are in genomic coords, not transcriptomic
        """
            
        reads_by_name = {}
        
        #for r in bamfile.fetch("chr1",24018269,24022915):
        #for r in bamfile.fetch("chr9",130209953,130213711):
        #for r in bamfile.fetch("chr10",99092254,99094458):
        for r in bamfile.fetch(contig):
            
            #only take unique reads - STAR -> 255 for uniq
            if r.mapping_quality!=255: 
                continue
            
            #only consider paired reads mapped on this contig
            if ((not r.is_paired) or (r.mate_is_unmapped) or (not r.next_reference_name==contig)):
                continue
            
            if not r.qname in reads_by_name:
                reads_by_name[r.qname] = {"name":r.qname}
             
            if r.is_read1:
                read = "read1"
            else:
                read = "read2"

            reads_by_name[r.qname]["%s_pos"%read] = r.pos
            reads_by_name[r.qname]["%s_l"%read] = r.query_alignment_length
            reads_by_name[r.qname]["%s_strand"%read] = not(r.is_reverse)

        strands = []
        lens = []
        poses = []
        names = []
        
        for r_name, r_info in reads_by_name.items():
            #print(r_info)
            #pdb.set_trace()
            if not "read1_strand" in r_info or not "read2_strand" in r_info:
                print(r_info)
                pdb.set_trace()
            
            strands.append(r_info['read1_strand'])
            names.append(r_name)
            
            assert r_info['read1_strand']!=r_info['read2_strand']

            if pe_3p:
                if r_info['read1_strand']:
                    poses.append(r_info['read2_pos']+r_info['read2_l'])
                    lens.append(-1)
                else:
                    poses.append(r_info['read2_pos'])
                    lens.append(-1)
            else:
                if r_info['read1_strand']:
                    poses.append(r_info['read1_pos'])
                    lens.append(-1)
                else:
                    poses.append(r_info['read1_pos']+r_info['read1_l'])
                    lens.append(-1)
        
        return strands, lens, poses, names


    def populate(self, bamfile, nosoftclipping, is_pe, pe_3p):
        """
        stores the position of the 5' and  3' end of each read
        for paired reads stores the 5' end and the 3' end of the other pair
        default is se
        """
        contigs = {x[0]:x[1] for x in zip(bamfile.references, bamfile.lengths)}
        
        F = open("%s.bed"%(self.fn),'w')
        for contig, clen in contigs.items():
            if "chrUn" in contig: continue
            if "random" in contig: continue
            
            #if contig!="chr10": continue ####TESTING

            sys.stderr.write("loading {contig}...".format(contig=contig))

            if is_pe:
                strands, lens, poses, names = self.load_pe(contig, clen, bamfile, pe_3p)
            else:
                strands, lens, poses = load_se(contig, clen, bamfile)


            l = len(poses)
            self.pos.append(np.array(poses,dtype='uint32'))
            self.strand.append(np.array(strands,dtype='int8'))
            self.length.append(np.array(lens,dtype='uint8'))
            self.contig_idx.append(np.ones(l)*self.curr_idx)
            
            self.contig_list.append(contig) 
            self.curr_idx+=1
            
            """ 
            for i, pos in enumerate(poses):
                l = strands[i] == 1 and 10 or -10
                s, e = sorted([pos, pos+l])
                F.write("{contig}\t{start}\t{end}\t{name}\n".format(contig=contig, start=s, end=e, name = names[i]))
            """ 
        
        F.close()

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
         
    def load_offsets(self, fn, no_offsets, read_length_range = None):
        
        if no_offsets:
            assert fn is None
            self.offsets = {}
            uniq_lens = np.unique(self.length)
            
            if read_length_range is not None:
                lo, hi = read_length_range
                uniq_lens = list(set(uniq_lens).intersection(set(range(lo, hi+1))))

            for l in uniq_lens:
                self.offsets[l] = 0
        else:
            t = pd.read_csv(fn, sep=" ", header=0)
            self.offsets = {}
            read_len_range = None
            if read_length_range is not None: 
                read_len_range = range(read_length_range[0], read_length_range[1]+1) 
            for i, row in t.iterrows():
                assert row['manual'] != -99
                if read_len_range is not None and not row['len'] in read_len_range:
                    continue
                self.offsets[row['len']] = row['manual']
        
    def get_offset_contig_arrays(self, contig, alignment_end):
        assert self.offsets is not None, "need to assign offsets"        
        
        contig_pos = self.pos[self.contig_idx == self.contig_to_idx[contig]]
        contig_strand = self.strand[self.contig_idx == self.contig_to_idx[contig]]
        contig_length = self.length[self.contig_idx == self.contig_to_idx[contig]]
        
        offset_contig_pos = np.zeros(contig_pos.shape[0],dtype='uint32')
        for l, o in self.offsets.items():
            
            alignment_end_offset = 0
            if alignment_end == "3p":
                alignment_end_offset = l 
            elif alignment_end == "5p":
                alignment_end_offset = 0
            else:
                assert False, "not supposed to happen!"

            w_pos = np.where((contig_length==l)&(contig_strand==1))
            w_neg = np.where((contig_length==l)&(contig_strand==0))
                    
            offset_contig_pos[w_pos] = contig_pos[w_pos]+o+alignment_end_offset
            offset_contig_pos[w_neg] = contig_pos[w_neg]-o-alignment_end_offset
        
        #ONLY KEEP THOSE THAT are amongst the lengths we want 
        #!!! need to assign strand before I mutate offset_contig_pos!
        contig_strand = contig_strand[offset_contig_pos!=0]
        offset_contig_pos = offset_contig_pos[offset_contig_pos!=0]
        
        s_neg = np.where(contig_strand==0)
        s_pos = np.where(contig_strand==1)
        
        return {0:offset_contig_pos[s_neg],
                1:offset_contig_pos[s_pos]}
    
    def get_offset_counts_by_contig(self, contig, alignment_end = "5p"):
        pos_arrays = self.get_offset_contig_arrays(contig, alignment_end = alignment_end)
        
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
    new_h5.populate(bamfile, args.nosoftclipping, args.is_pe, args.pe_3p)
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



def get_coverage_info(counts_by_contig, cvg_ob):
    
    if cvg_ob.strand == 1:
        start, stop = cvg_ob.coding_exons[0][0], cvg_ob.coding_exons[-1][1] 
        start_e, stop_e = [[start,start+3]],[[stop-3, stop]]
    else:
        start, stop = cvg_ob.coding_exons[-1][1]-1, cvg_ob.coding_exons[0][0]
        start_e, stop_e = [[start-2,start+1]], [[stop, stop+3]]
    
    strand = cvg_ob.strand 
    CDS_cvg, CDS_cvg_u, CDS_l = get_cvg(cvg_ob.coding_exons, strand, counts_by_contig)
    CDS_l200_cvg, CDS_l200_cvg_U, CDS_l200_l = get_cvg(cvg_ob.coding_exons, strand, counts_by_contig, maxlen=-200)
    
    UTR_3p_cvg, UTR_3p_cvg_u, UTR_3p_l = get_cvg(cvg_ob.UTR_3p_exons, strand, counts_by_contig)
    UTR_3p200_cvg, UTR_3p200_cvg_U, UTR_3p200_l = get_cvg(cvg_ob.UTR_3p_exons, strand, counts_by_contig, maxlen=200)

    UTR_5p_cvg, UTR_5p_cvg_u, UTR_5p_l = get_cvg(cvg_ob.UTR_5p_exons, strand, counts_by_contig)
    STOP_cvg, CDS_cvg_U, l = get_cvg(stop_e, strand, counts_by_contig)
    START_cvg, CDS_cvg_U, l = get_cvg(start_e, strand, counts_by_contig)
     
    return {"CDS_cvg":CDS_cvg,
            "CDS_len":CDS_l,
            "CDS_up200_cvg":CDS_l200_cvg,
            "CDS_up200_l":CDS_l200_l,
            "UTR_5p_cvg":UTR_5p_cvg,
            "UTR_5p_len":UTR_5p_l,
            "UTR_3p_cvg":UTR_3p_cvg,
            "UTR_3p_len":UTR_3p_l,
            "UTR_3p200_cvg":UTR_3p200_cvg,
            "UTR_3p200_l":UTR_3p200_l,
            "START_cvg":START_cvg,
            "STOP_cvg":STOP_cvg}

def get_cvg(exons, strand, counts_by_contig, maxlen=None):
    l = 0
    cvg = 0
    u_cvg = 0
    cvg_array = []
    u_cvg_array = []
    for e in exons:
        s,e = e
        l+=e-s
        w = np.where((counts_by_contig[strand]['pos']>=s)&(counts_by_contig[strand]['pos']<e))
        
        cvg_array.append(counts_by_contig[strand]['count'][w])
        u_cvg_array.append(counts_by_contig[strand]['count'][w]!=0)
        
        cvg += np.sum(cvg_array[-1])
        u_cvg = np.sum(u_cvg_array[-1])
    
    if maxlen is not None and cvg>0:
        cvg_array = np.concatenate(cvg_array)
        u_cvg_array = np.concatenate(u_cvg_array)
        if abs(maxlen) <l:
            l = abs(maxlen)
            if maxlen>0:
                cvg = np.sum(cvg_array[:maxlen])
                u_cvg = np.sum(u_cvg_array[:maxlen])
            else:
                cvg = np.sum(cvg_array[maxlen:])
                u_cvg = np.sum(u_cvg_array[maxlen:])
    return cvg, u_cvg, l

    """
    ACTUALLY SLOWEWR!!!
    delta_t_1 = timeit.default_timer() - t1
    if len(exons) > 0:
        t1 = timeit.default_timer()
        poses = np.concatenate([np.arange(e[0], e[1]) for e in exons]) 
        bool_poses = np.in1d(counts_by_contig[strand]['pos'], poses)
        cvg_2 = np.sum(counts_by_contig[strand]['count'][bool_poses])
        delta_t_2 = timeit.default_timer() - t1
        print(cvg, cvg_2, delta_t_1/delta_t_2)
    """
            

###############
def get_binned_cvg(exons, strand, counts_by_contig):
    l = 0
    cvgs = []
    for e in exons:
        s,e = e
        l+=e-s
        w = np.where((counts_by_contig[strand]['pos']>=s)&(counts_by_contig[strand]['pos']<e))
        cvgs.append(counts_by_contig[strand]['count'][w])
    return cvg, l

def get_bp_coverage(counts_by_contig, cvg_ob, width, feature_pos):
    
    strand = cvg_ob.strand
    s,e = feature_pos-width, feature_pos+width+1
    cvg_vect = np.zeros(2*width+1)
    u_cvg_vect = np.zeros(2*width+1)
    
    w = np.where((counts_by_contig[strand]['pos']>=s)&(counts_by_contig[strand]['pos']<e))
    poses = counts_by_contig[strand]['pos'][w]
    cvg = counts_by_contig[strand]['count'][w]
    u_cvg = (cvg!=0)*1

    if strand==1: 
        centered_poses = poses-feature_pos
    else:
        centered_poses = feature_pos-poses
    
    cvg_vect[centered_poses+width] = cvg
    u_cvg_vect[centered_poses+width] = u_cvg
    
    return cvg_vect, u_cvg_vect
    
    
def RFPCountTable(args):
    h5 = h5_ribo(args.fn_h5)
    h5.load_offsets(args.fn_aSiteOffsets, 
                    args.no_aSiteOffsets, 
                    args.read_length_range)
    
    gene_id_subset = None 
    if args.fn_gene_subset is not None:
        t_subset = pd.read_csv(args.fn_gene_subset, header=0, sep="\t")
        gene_id_subset = t_subset.gene_id.values
    
    feature_table = None 
    if args.feature=="OTHER":
        feature_table = pd.read_csv(args.fn_features, header=0, sep="\t", index_col=["transcript_id","feature_idx"])

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
        t = timeit.default_timer()
        counts_by_contig = h5.get_offset_counts_by_contig(contig, alignment_end = args.alignment_end)
        
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
            elif args.feature=="OTHER":
                if not (cvg_ob.TID,0) in feature_table.index:
                    continue
                pos = feature_table.ix[(cvg_ob.TID,0)]['g_start']
            else:
                assert False, "no method for feature: %s"%(args.feature) 

            cvg_vect, u_cvg_vect =  get_bp_coverage(counts_by_contig, cvg_ob, args.width, pos)
            cvg_info = get_coverage_info(counts_by_contig, cvg_ob)
                
            #######
            #HERE add the total sum, and 1,2,3 counts
            cum_prop = [0,0,0,0,0]
            if np.sum(cvg_vect)>0:
                k_largest_inds = np.argpartition(cvg_vect,-5)[-5:]
                k_largest_vals = np.sort(cvg_vect[k_largest_inds])[::-1]
                s = np.sum(cvg_vect) 
                cum_prop = np.cumsum(k_largest_vals)/s
                
            for i, pos in enumerate(range(-args.width,args.width)):
                outrows.append({"tid":cvg_ob.TID,
                                "gene_id":cvg_ob.gene_id,
                                "gene_name":cvg_ob.g.names[0],
                                "cvg":cvg_vect[i],
                                "u_cvg":u_cvg_vect[i],
                                "CDS_cvg":cvg_info['CDS_cvg'],
                                "CDS_len":cvg_info['CDS_len'],
                                "pos":pos,
                                "feature":args.feature,
                                "cumprop_1":cum_prop[0],
                                "cumprop_2":cum_prop[1]})

        print("1", timeit.default_timer()-t, len(cvg_objs))
    print("done")
    t_time = timeit.default_timer()
    t = pd.DataFrame(outrows)
    print("time to make dataframe:", timeit.default_timer()-t_time)
    #t.to_csv(args.fn_out, sep="\t", index=False, compression="gzip")
    t.to_hdf(args.fn_out, key="data", mode="w", format='t', complib='zlib', complevel=5)
    print("time to output hdf:", timeit.default_timer()-t_time)


def makeFeatureSummary(args):
    if args.fn_features is not None:
        makeFeatureSummaryFromTable(args)
    else:
        makeGeneFeatureSummary(args)

def makeFeatureSummaryFromTable(args):
    """
    get coverage over all the features in a table
    each row must include:
        contig
        start
        end
        strand
    maintain all other columns in the table
    """
    h5 = h5_ribo(args.fn_h5)
    h5.load_offsets(args.fn_aSiteOffsets, 
                    args.no_aSiteOffsets)

    t_features = pd.read_csv(args.fn_features, header=0, sep="\t")
    cvgs = []
    
    s_key, e_key = args.feature_start_key, args.feature_end_key

    for contig, features in t_features.groupby("contig"):
        if not contig in h5.contig_to_idx: continue
        
        sys.stderr.write("{contig}...".format(contig=contig))
        counts_by_contig = h5.get_offset_counts_by_contig(contig) #pos #count
        
        for i, feature_row in features.iterrows():
            s,e,strand = (feature_row[s_key], 
                          feature_row[e_key],
                          feature_row['strand'])

            cvg, u_cvg, l = get_cvg([[s,e]], strand, counts_by_contig)
            ###########
            cvg_vect, u_cvg_vect =  get_bp_coverage(counts_by_contig, cvg_ob, args.width, pos)
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

            #######
            cvgs.append({"contig":contig, 
                         s_key:s,
                         e_key:e,
                         "cvg":cvg})
    t_cvg = pd.DataFrame(cvgs)    
    T = pd.merge(t_features, 
                 t_cvg, 
                 left_on=['contig', s_key, e_key], 
                 right_on = ['contig', s_key, e_key])
    
    T.to_csv(args.fn_out, sep="\t", index=False, compression="gzip")
    
def makeGeneFeatureSummary(args):

    h5 = h5_ribo(args.fn_h5)
    h5.load_offsets(args.fn_aSiteOffsets, 
                    args.no_aSiteOffsets)

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
            
            for i, e in enumerate(cvg_ob.coding_exons):
                cvg, u_cvg, l = get_cvg([e], cvg_ob.strand, counts_by_contig)
                d = {"type":"exon",
                     "type_idx":i,
                     "cvg":cvg,
                     "len":l}
                d.update(cvg_inf)
                outrows.append(d)
                
            for i, intr in enumerate(get_introns(cvg_ob.coding_exons)):
                cvg, u_cvg, l = get_cvg([intr], cvg_ob.strand, counts_by_contig)
                d = {"type":"intron",
                     "type_idx":i,
                     "cvg":cvg,
                     "len":l}
                d.update(cvg_inf)
                outrows.append(d)

    t = pd.DataFrame(outrows)
    t.to_csv(args.fn_out, sep="\t", index=False, compression="gzip")

def makeSummary(args):

    h5 = h5_ribo(args.fn_h5)
    h5.load_offsets(args.fn_aSiteOffsets, 
                    args.no_aSiteOffsets, 
                    args.read_length_range)
    
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


def get_flat_gene_regions(cvg_objs_by_contig):
    """
    only protein coding
    """
    flat_regions_by_contig = {}
    for contig, cvg_objs in cvg_objs_by_contig.items():
        
        sys.stderr.write("{contig}...".format(contig=contig))
        arrays = []
        
        for cvg_ob in cvg_objs:
            arrays.append(np.concatenate([np.arange(e[0],e[1]) for e in cvg_ob.exons]))

        flat_regions_by_contig[contig] = np.unique(np.concatenate(arrays))
    
    return flat_regions_by_contig

def makeReadLengthSummary(args):

    h5 = h5_ribo(args.fn_h5)
    h5.load_offsets(args.fn_aSiteOffsets, args.no_aSiteOffsets)
    
    logger = logging.getLogger(args.fn_logfile)
    s_id, path, genes = get_indexed_genes_for_identifier(args.fn_gtf_index,
                                                         logger, 
                                                         args.gtf_ID)
    sys.stderr.write("done\n")
    sys.stderr.write("loading gtf\n")
    cvg_objs_by_contig = get_cvg_objs_by_contig(genes,
                                                "transcript")
    sys.stderr.write("done\n")
    
    flat_regions_by_contig = get_flat_gene_regions(cvg_objs_by_contig)
     
    coding_counts_by_size = {} 
    for contig, flat_regions in flat_regions_by_contig.items():
        if not contig in h5.contig_to_idx: continue

        curr_idx = h5.contig_to_idx[contig]
        w_curr_contig = np.where(h5.contig_idx==curr_idx)
        poses = h5.pos[w_curr_contig]
        lengths = h5.length[w_curr_contig]
        
        bool_coding_regions = np.in1d(poses, flat_regions)
        coding_lens = lengths[bool_coding_regions] 
        
        U, counts = np.unique(coding_lens, return_counts = True)
        for i in range(U.shape[0]):
            l = U[i]
            count = counts[i]
            if not l in coding_counts_by_size:
                coding_counts_by_size[l] = 0
            coding_counts_by_size[l]+=count
    
    """
    simple summary inf
    """

    U, counts = np.unique(h5.length, return_counts = True)
    outrows = []
    for i in range(U.shape[0]):
        l = U[i]
        count = counts[i]
        ####
        U_pos = np.unique(h5.pos[h5.length==l])
        n_unique_pos = U_pos.shape[0] 
        
        if not l in coding_counts_by_size:
            coding_counts_by_size[l] = 0

        outrows.append({"length":l,
                        "count":count,
                        "coding_counts":coding_counts_by_size[l],
                        "unique_positions":n_unique_pos})
        
    t = pd.DataFrame(outrows)
    t.to_csv(args.fn_out, sep="\t", index=False)


def makeWig(args):
    
    h5 = h5_ribo(args.fn_h5)
    h5.load_offsets(args.fn_aSiteOffsets, 
                    args.no_aSiteOffsets)
    
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
    parser_create.add_argument("--nosoftclipping", action='store_true', default=False)
    parser_create.add_argument("--is_pe", action='store_true', default=False)
    parser_create.add_argument("--pe_3p", action='store_true', default=False)
    parser_create.set_defaults(func=create)
    
    #make wig from h5
    parser_makeWig = subparsers.add_parser("makeWig")
    parser_makeWig.add_argument("--fn_h5", required=True)
    parser_makeWig.add_argument("--fn_out", required=True)
    parser_makeWig.add_argument("--fn_aSiteOffsets", required=True)
    parser_makeWig.add_argument("--no_aSiteOffsets", required=False, 
                                                     default=False, 
                                                     action='store_true')
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
    parser_makeSum.add_argument("--fn_aSiteOffsets", required=False)
    parser_makeSum.add_argument("--no_aSiteOffsets", required=False, 
                                                     default=False, 
                                                     action='store_true')
    parser_makeSum.add_argument("--read_length_range", 
                                                     required=False, 
                                                     default=None,
                                                     type=int,
                                                     nargs=2)
    parser_makeSum.add_argument("--fn_gtf_index", required=True)
    parser_makeSum.add_argument("--gtf_ID", required=True)
    parser_makeSum.add_argument("--fn_logfile", default='/dev/stderr')
    parser_makeSum.set_defaults(func=makeSummary)
    
    #make summary of readlen details
    parser_makeSum = subparsers.add_parser("makeReadLengthSummary")
    parser_makeSum.add_argument("--fn_h5", required=True)
    parser_makeSum.add_argument("--fn_out", required=True)
    parser_makeSum.add_argument("--fn_aSiteOffsets", required=True)
    parser_makeSum.add_argument("--no_aSiteOffsets", required=False, 
                                                     default=False, 
                                                     action='store_true')
    parser_makeSum.add_argument("--fn_gtf_index", required=True)
    parser_makeSum.add_argument("--gtf_ID", required=True)
    parser_makeSum.add_argument("--fn_logfile", default='/dev/stderr')
    parser_makeSum.set_defaults(func=makeReadLengthSummary)
    
    
    #make summary of coverage details by feature
    parser_makeFeatureSum = subparsers.add_parser("makeFeatureSummary")
    parser_makeFeatureSum.add_argument("--fn_h5", required=True)
    parser_makeFeatureSum.add_argument("--fn_out", required=True)
    parser_makeFeatureSum.add_argument("--fn_aSiteOffsets", required=False)
    parser_makeFeatureSum.add_argument("--no_aSiteOffsets", required=False, 
                                                     default=False, 
                                                     action='store_true')
    parser_makeFeatureSum.add_argument("--fn_features", required=False)
    parser_makeFeatureSum.add_argument("--feature_start_key", 
                                       required=False, 
                                       default="start")
    parser_makeFeatureSum.add_argument("--feature_end_key", 
                                       required=False, 
                                       default="end")
    parser_makeFeatureSum.add_argument("--fn_gtf_index", required=False)
    parser_makeFeatureSum.add_argument("--gtf_ID", required=False)
    parser_makeFeatureSum.add_argument("--fn_logfile", default='/dev/stderr')
    parser_makeFeatureSum.set_defaults(func=makeFeatureSummary)
    
    #output feature centered counts
    parser_makeRFPCountTable = subparsers.add_parser("RFPCountTable")
    parser_makeRFPCountTable.add_argument("--fn_h5", required=True)
    parser_makeRFPCountTable.add_argument("--fn_out", required=True)
    parser_makeRFPCountTable.add_argument("--fn_aSiteOffsets", required=False)
    parser_makeRFPCountTable.add_argument("--no_aSiteOffsets", required=False, 
                                                     default=False, 
                                                     action='store_true')
    parser_makeRFPCountTable.add_argument("--read_length_range", 
                                                     required=False, 
                                                     default=None,
                                                     type=int,
                                                     nargs=2)
    parser_makeRFPCountTable.add_argument("--fn_gtf_index", required=True)
    parser_makeRFPCountTable.add_argument("--gtf_ID", required=True)
    parser_makeRFPCountTable.add_argument("--fn_logfile", default='/dev/stderr')
    parser_makeRFPCountTable.add_argument("--feature", required=True, choices = ["START","STOP","POLYA","OTHER"])
    parser_makeRFPCountTable.add_argument("--fn_features", required=False, default=None)
    parser_makeRFPCountTable.add_argument("--alignment_end", required=True, choices = ["5p","3p"], default='5p')
    parser_makeRFPCountTable.add_argument("--width", default=50, type=int)
    parser_makeRFPCountTable.add_argument("--fn_gene_subset", required=False, default=None)
    parser_makeRFPCountTable.set_defaults(func=RFPCountTable)
    
    args = parser.parse_args()
    args.func(args)


