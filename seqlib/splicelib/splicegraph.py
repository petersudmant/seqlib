import numpy as np
import pdb
import argparse
import BCBio.GFF as GFF
import re
import string 
import networkx as nx

from transcript import Transcript
from collections import defaultdict
from sys import stderr
from fastahack import FastaHack

from bx.intervals.intersection import Interval, IntervalTree

trans = string.maketrans('ATCGatcg', 'TAGCtacg')

def revcomp(s):
    return s.translate(trans)[::-1]


class SpliceGraph(object):

    def __init__(self, contig, **kwargs):
            
        """
          5'      3'
        + GT------AG
        
          3'      5'
        - CT------AC

        a collection of hashes of lists representing all splice site 
        and all exon boundary relationships

        """
        
        self.contig = contig
        
        """
        exon-intron-exon junctions
        """
        self.F_3p_5p_ss = kwargs["F_3p_5p_ss"]
        self.F_5p_3p_ss = kwargs["F_5p_3p_ss"]
        self.R_3p_5p_ss = kwargs["R_3p_5p_ss"]
        self.R_5p_3p_ss = kwargs["R_5p_3p_ss"]
        
        """
        intron-exon-intron junctions
        """
        self.F_exon_s_e = kwargs["F_exon_s_e"]
        self.F_exon_e_s = kwargs["F_exon_e_s"]
        self.R_exon_s_e = kwargs["R_exon_s_e"]
        self.R_exon_e_s = kwargs["R_exon_e_s"]
       
        """
        last exon s-e/e-s junctions
        """
        self.F_le_s_e = kwargs["F_le_s_e"]
        self.F_le_e_s = kwargs["F_le_e_s"]
        self.R_le_s_e = kwargs["R_le_s_e"]
        self.R_le_e_s = kwargs["R_le_e_s"]
        
        """
        first exon s-e/e-s junctions
        """
        self.F_fe_s_e = kwargs["F_fe_s_e"]
        self.F_fe_e_s = kwargs["F_fe_e_s"]
        self.R_fe_s_e = kwargs["R_fe_s_e"]
        self.R_fe_e_s = kwargs["R_fe_e_s"]

        """
        5'/3' ss to gene_info
        """
        self.F_5p_to_gene_info = kwargs["F_5p_to_gene_info"]
        self.F_3p_to_gene_info = kwargs["F_3p_to_gene_info"]
        self.R_5p_to_gene_info = kwargs["R_5p_to_gene_info"]
        self.R_3p_to_gene_info = kwargs["R_3p_to_gene_info"]
        
        """
        interval tree of introns with introns represented from "left" to "right"
        regardless of strand, ie, s<e
        """
        self.F_LR_intron_interval_tree = kwargs["F_LR_intron_interval_tree"]
        self.R_LR_intron_interval_tree = kwargs["R_LR_intron_interval_tree"]
        
        """
        exon graphs for easy connected compenent finding
        all exons stored in s<e format irrespective of strand
        """
        self.F_exon_G = kwargs["F_exon_G"]
        self.R_exon_G = kwargs["R_exon_G"]
        
        """
        nx based ss graph with 
        """
        self.ss_G = kwargs["ss_G"]
        
        self.all_uniq_exons = kwargs["all_uniq_exons"]
        
    def get_csx(self, ss, strand, ss_type_3p=False, ss_type_5p=False):
        return self.get_common_shortest_exon(ss, 
                                             strand, 
                                             ss_type_3p=ss_type_3p, 
                                             ss_type_5p=ss_type_5p)

    def get_common_shortest_exon(self, ss, strand, ss_type_3p=False, ss_type_5p=False):
        """
        ss_type = 5' | 3'
        strand = FWD | REV 

        all exons should be returned s<e

        NOTE: given 
                5'-1
            -----  5'-2
            --------
            this will return 5'-1
            BUT, given
            
                5'-1
            -----  5'-2
            --------
              -- 5'-3
            this will still return 5'-1, I THINK this is correct...
        """
        assert ss_type_3p!=ss_type_5p, "must pick 3' or 5' (_3p, _5p)"
        assert strand=="FWD" or strand=="REV", "must pick FWD or REV (FWD, REV)"

        if strand=="FWD":
            if ss_type_5p:
                return tuple([max(self.F_exon_e_s[ss]), ss])
            else:
                return tuple([ss, min(self.F_exon_s_e[ss])])
        else:
            if ss_type_5p:
                return tuple([ss, min(self.R_exon_s_e[ss])])
            else:
                return tuple([max(self.R_exon_e_s[ss]), ss])
    
    def get_first_last_exons(self, strand):
        if strand=="REV":
            le_s_e = self.R_le_s_e
            le_e_s = self.R_le_e_s
            fe_s_e = self.R_fe_s_e
            fe_e_s = self.R_fe_e_s
        else:
            le_s_e = self.F_le_s_e
            le_e_s = self.F_le_e_s
            fe_s_e = self.F_fe_s_e
            fe_e_s = self.F_fe_e_s

        return le_s_e, le_e_s, fe_s_e, fe_e_s

    def get_intron_exon_juncs(self, strand):
        if strand=="REV":
            i_5p_3p = self.R_5p_3p_ss
            i_3p_5p = self.R_3p_5p_ss
            e_3p_5p = self.R_exon_e_s
            e_5p_3p = self.R_exon_s_e
        else:
            i_5p_3p = self.F_5p_3p_ss
            i_3p_5p = self.F_3p_5p_ss
            e_3p_5p = self.F_exon_s_e
            e_5p_3p = self.F_exon_e_s

        return i_5p_3p, i_3p_5p, e_5p_3p, e_3p_5p
    
    def get_gene_ID(self, strand_d, ss_3p=None, ss_5p=None):
        
        assert ss_3p!=ss_5p

        if strand_d == "FWD":
            if ss_3p:
                return self.F_3p_to_gene_info[ss_3p][0]['gene_ID']
            else:
                return self.F_5p_to_gene_info[ss_5p][0]['gene_ID']
        else:
            if ss_3p:
                return self.R_3p_to_gene_info[ss_3p][0]['gene_ID']
            else:
                return self.R_5p_to_gene_info[ss_5p][0]['gene_ID']
    
    def get_gene_name(self, strand_d, ss_3p=None, ss_5p=None):
        
        assert ss_3p!=ss_5p
        try:
            if strand_d == "FWD":
                if ss_3p:
                    return self.F_3p_to_gene_info[ss_3p][0]['gene_name']
                else:
                    return self.F_5p_to_gene_info[ss_5p][0]['gene_name']
            else:
                if ss_3p:
                    return self.R_3p_to_gene_info[ss_3p][0]['gene_name']
                else:
                    return self.R_5p_to_gene_info[ss_5p][0]['gene_name']
        except:
            return None
                
    def get_gene_info(self, strand_d, ss_3p, ss_5p):

        if strand_d == "FWD":
            us_gene_inf = self.F_5p_to_gene_info[ss_5p]
            ds_gene_inf = self.F_3p_to_gene_info[ss_3p]
        else:
            us_gene_inf = self.R_5p_to_gene_info[ss_5p]
            ds_gene_inf = self.R_3p_to_gene_info[ss_3p]

        gene_inf = []
        for inf in us_gene_inf+ds_gene_inf:
            if not inf in gene_inf: 
                gene_inf.append(inf)
        return gene_inf

    def get_skipped_alt_exon(self, seq, F_gff, F_novel_gff, F_bed, junc_writer, n=60):
        """
        get all 5' 3' ss pairs and search in-between for lil' guys
        """

        for strand_d, ss_juncs in { "FWD" : self.F_5p_3p_ss, "REV" : self.R_5p_3p_ss }.iteritems():
            t_complete=0 
            for annot_ss_5p, annot_ss_3p_list in ss_juncs.iteritems():
                
                discovered_alt_transcript_hashes = {}
                
                t_complete+=1
                if t_complete%100 ==0 :
                    print t_complete, "/", len(ss_juncs.keys())

                for annot_ss_3p in annot_ss_3p_list: 
                    
                    ds_exon = self.get_common_shortest_exon(annot_ss_3p, strand_d, ss_type_3p=True)
                    us_exon = self.get_common_shortest_exon(annot_ss_5p, strand_d, ss_type_5p=True)
                
                    """
                    GET GENE INFO
                    """
                    gene_inf = self.get_gene_info(strand_d, annot_ss_3p, annot_ss_5p)
                
                    if strand_d=="FWD":
                        ss_donor_seq="GT"
                        ss_acceptor_seq="AG"
                        delta_donor=0
                        delta_acceptor=2
                        seq_s, seq_e = annot_ss_5p+2, annot_ss_3p-2
                    else:
                        ss_donor_seq="AC"
                        ss_acceptor_seq="CT"
                        delta_donor=2
                        delta_acceptor=0
                        seq_s, seq_e = annot_ss_3p+2, annot_ss_5p-2 
                    

                    donors = np.array([m.start()+seq_s+delta_donor for m in re.finditer(ss_donor_seq, seq[seq_s:seq_e].upper())])
                    acceptors = np.array([m.start()+seq_s+delta_acceptor for m in re.finditer(ss_acceptor_seq, seq[seq_s:seq_e].upper())])
                    
                    """
                    mult lets us just get the exons in the right direction
                    """
                    
                    if strand_d=="FWD": 
                        mult=1
                    else:
                        mult=-1

                    for donor in donors:
                        ds = mult*(donor-acceptors)
                        for acceptor in acceptors[(ds<n)&(ds>0)]:
                            source = ""
                            alt_exon=sorted([acceptor, donor])
                        
                            if (alt_exon[1]==alt_exon[0]):
                                continue
                            exon_paths = {"A":[0,2], "B":[0,1,2]}

                            EXONS  = [us_exon, alt_exon, ds_exon]
                            """
                            FEATURE_ID = "%s_%s"%(source, "_".join(["%s:%d-%d"%(self.contig, e[0], e[1]) for e in EXONS]))
                            GENE_NAME = "%s_%s"%(source, ",".join(["%s"%gi['gene_name'] for gi in gene_inf]))
                            GENE_ID = "%s_%s"%(source,"_".join(["%s"%gi['gene_ID'] for gi in gene_inf]))
                            """
                            FEATURE_ID, GENE_NAME, GENE_ID="","",""
                            G_START = min(us_exon[0], ds_exon[0])
                            G_END = max(us_exon[1], ds_exon[1])
                            STRAND = strand_d == "FWD" and 1 or -1
                        
                            alt_T = Transcript(contig = self.contig, 
                                               feature_ID = FEATURE_ID,
                                               exons = EXONS, 
                                               gene_name = GENE_NAME,
                                               gene_ID = GENE_ID,
                                               g_start = G_START,
                                               g_end = G_END,
                                               strand = STRAND)

                            #T_hash = alt_T.get_tuple_hash()
                            #discovered_alt_transcript_hashes[T_hash] = 1

                            junc_tups = alt_T.junc_tuples(exon_paths, source, True)
                            junc_writer.write(junc_tups)
                            """ 
                            gff_s = alt_T.gff_string(exon_paths, source)
                            gff_novel_s = alt_T.gff_string({"A":[0,1]}, source)
                            bed_s = alt_T.bed_string(exon_paths, source, True)
                            
                            F_gff.write(gff_s)
                            F_novel_gff.write(gff_novel_s)
                            F_bed.write(bed_s)
                            """ 


    def get_5p_micro_exon(self, seq, F_gff, F_novel_gff, F_bed, junc_writer, micro_exon_dir, n=60):
        get_5p_3p_alt_exon(self, seq, F_gff, F_novel_gff, F_bed, junc_writer, get_5p=True, n=n)
    
    def get_3p_micro_exon(self, seq, F_gff, F_novel_gff, F_bed, junc_writer, micro_exon_dir, n=60):
        get_5p_3p_alt_exon(self, seq, F_gff, F_novel_gff, F_bed, junc_writer, get_3p=True, n=n)

    def get_5p_3p_alt_exon(self, seq, F_gff, F_novel_gff, F_bed, junc_writer, get_5p=False, get_3p=False, n=60):
        """
        get all putative alternative 5'/3' alternative exons 
        """
        assert get_5p != get_3p, "EITHER get_5p or get_3p (5' or 3') must be passed"

        if get_3p:
            fwd_ss_juncs = self.F_5p_3p_ss
            rev_ss_juncs = self.R_5p_3p_ss
        else:
            fwd_ss_juncs = self.F_3p_5p_ss
            rev_ss_juncs = self.R_3p_5p_ss

        for strand_d, ss_juncs in { "FWD" : fwd_ss_juncs, "REV" : rev_ss_juncs }.iteritems():
            """
            whether we are looking for alt 5' or 3' ss, one will be held constant while 
            we scan around the other, call the one which we hold constant to be "fixed"
            and the other to be called alt
            """
            for fixed_ss_3_or_5, alt_annot_ss_5_or_3_list in ss_juncs.iteritems():
                
                discovered_alt_transcript_hashes = {}

                for alt_annot_5_or_3 in alt_annot_ss_5_or_3_list: 
                    
                    if get_3p:
                        curr_3p = alt_annot_5_or_3
                        curr_5p = fixed_ss_3_or_5
                    if get_5p:
                        curr_3p = fixed_ss_3_or_5   
                        curr_5p = alt_annot_5_or_3

                    ds_exon = self.get_common_shortest_exon(curr_3p, strand_d, ss_type_3p=True)
                    us_exon = self.get_common_shortest_exon(curr_5p, strand_d, ss_type_5p=True)
                
                    gene_inf = self.get_gene_info(strand_d, curr_3p,curr_5p)

                    if get_3p:
                        ss_seq = strand_d == "FWD" and "AG" or "CT"
                    else:
                        ss_seq = strand_d == "FWD" and "GT" or "AC"
                            
                    min_alt, max_alt = alt_annot_5_or_3-n, alt_annot_5_or_3+n
                    """
                    delta offsets are to adjust the identified splice donor/acceptor to 
                    the start/end
                    """
                    if strand_d == "FWD":
                        if get_5p:
                            seq_s, seq_e = max(min_alt, us_exon[0]), min(max_alt, ds_exon[0]-2)
                            delta=0
                        else: 
                            seq_s, seq_e = max(min_alt, us_exon[1]+2), min(max_alt, ds_exon[1])
                            delta=2
                    else:
                        if get_5p:
                            seq_s, seq_e = max(min_alt, ds_exon[1]+2), min(max_alt, us_exon[1])
                            delta=2
                        else:
                            seq_s, seq_e = max(min_alt, ds_exon[0]), min(max_alt, us_exon[0]-2)
                            delta=0
                    """
                    THIS IS RUNNING MULTIPLE TIMES MORE THANT IT SHOULD REALLY... 
                    """
                    alt_ss_list = np.array([m.start()+seq_s+delta for m in re.finditer(ss_seq, seq[seq_s:seq_e].upper())])

                    for alt_ss in alt_ss_list:
                        if alt_ss == alt_annot_5_or_3: 
                            continue

                        if alt_ss in alt_annot_ss_5_or_3_list: 
                            source = "annotated"
                        else:
                            source = "novel"
                        
                        if get_3p:
                            exon_paths = {"A":[0,1], "B":[0,2]}
                        else:
                            exon_paths = {"A":[0,2], "B":[1,2]}

                        if strand_d=="FWD" and get_3p:
                            alt_exon=[alt_ss, ds_exon[1]]
                        elif strand_d=="REV" and get_3p:
                            alt_exon=[ds_exon[0], alt_ss]
                        elif strand_d=="FWD" and get_5p:
                            alt_exon=[us_exon[0], alt_ss]
                        elif strand_d=="REV" and get_5p:
                            alt_exon=[alt_ss, us_exon[1]]
                        
                        #no, 0 length exons
                        if alt_exon[1]==alt_exon[0]: continue

                        EXONS  = [us_exon, alt_exon, ds_exon]
                        FEATURE_ID = "%s_%s"%(source, "_".join(["%s:%d-%d"%(self.contig, e[0], e[1]) for e in EXONS]))
                        GENE_NAME = "%s_%s"%(source, ",".join(["%s"%gi['gene_name'] for gi in gene_inf]))
                        GENE_ID = "%s_%s"%(source,"_".join(["%s"%gi['gene_ID'] for gi in gene_inf]))
                        G_START = min(us_exon[0], ds_exon[0])
                        G_END = max(us_exon[1], ds_exon[1])
                        STRAND = strand_d == "FWD" and 1 or -1
                        
                        alt_T = Transcript(contig = self.contig, 
                                           feature_ID = FEATURE_ID,
                                           exons = EXONS, 
                                           gene_name = GENE_NAME,
                                           gene_ID = GENE_ID,
                                           g_start = G_START,
                                           g_end = G_END,
                                           strand = STRAND)

                        T_hash = alt_T.get_tuple_hash()
                        if not T_hash in discovered_alt_transcript_hashes:
                            discovered_alt_transcript_hashes[T_hash] = 1

                            gff_s = alt_T.gff_string(exon_paths, source)
                            gff_novel_s = alt_T.gff_string({"A":[0,1]}, source)
                            bed_s = alt_T.bed_string(exon_paths, source, True)
                            junc_tups = alt_T.junc_tuples(exon_paths, source, True)
                            
                            F_gff.write(gff_s)
                            F_novel_gff.write(gff_novel_s)
                            F_bed.write(bed_s)
                            junc_writer.write(junc_tups)

    def enumerate_splice_junctions(self, seq):
        fwd_5p = defaultdict(int)
        fwd_3p = defaultdict(int)
        rev_5p = defaultdict(int)
        rev_3p = defaultdict(int)

        for F_5p, F_3p_list in self.F_5p_3p_ss.iteritems():
            for F_3p in F_3p_list: 
                fwd_3p[seq[F_3p-2:F_3p]]+=1
        for F_3p, F_5p_list in self.F_3p_5p_ss.iteritems():
            for F_5p in F_5p_list: 
                fwd_5p[seq[F_5p:F_5p+2]]+=1
        for R_5p, R_3p_list in self.R_5p_3p_ss.iteritems():
            for R_3p in R_3p_list: 
                rev_3p[revcomp(seq[R_3p:R_3p+2])]+=1
        for R_3p, R_5p_list in self.R_3p_5p_ss.iteritems():
            for R_5p in R_5p_list: 
                rev_5p[revcomp(seq[R_5p-2:R_5p])]+=1 
        
        print "fwd"
        print fwd_5p, np.sum(fwd_5p.values())
        print fwd_3p, np.sum(fwd_3p.values())
        print "rev"
        print rev_5p, np.sum(rev_5p.values())
        print rev_3p, np.sum(rev_3p.values())


def init_splice_graphs_from_gff3(fn_gff, **kwargs):
    """
    limit_info={"gff_id": ["13"]}
    """
    contigs = kwargs.get('contigs', [])
    features = kwargs.get('features', ["transcript", "mRNA", "protein_coding"])
    min_exons = kwargs.get('min_exons', 1) 

    limit_info = {"gff_id": contigs}
    GFF_parser = GFF.GFFParser()
    
    SGs_by_contig = {}
    for rec in GFF_parser.parse_in_parts(open(fn_gff), limit_info=limit_info):
        #iterate over chromosomes
        contig = rec.id
        print >>stderr, "parsing records for record id:%s..."%rec.id
        #contig_seq = fa.get_sequence(rec.id)
                   
        F_3p_5p_ss = {} 
        F_5p_3p_ss = {}  
        R_3p_5p_ss = {} 
        R_5p_3p_ss = {} 

        F_exon_s_e = {} 
        F_exon_e_s = {} 
        R_exon_s_e = {} 
        R_exon_e_s = {} 
        
        F_le_s_e = {}
        F_le_e_s = {}
        R_le_s_e = {}
        R_le_e_s = {}
        
        F_fe_s_e = {}
        F_fe_e_s = {}
        R_fe_s_e = {}
        R_fe_e_s = {}

        F_LR_intron_interval_tree = IntervalTree()
        R_LR_intron_interval_tree = IntervalTree()
        F_added_LR_intron_intervals = {}
        R_added_LR_intron_intervals = {}

        F_5p_to_gene_info = {}
        F_3p_to_gene_info = {}
        R_5p_to_gene_info = {}
        R_3p_to_gene_info = {}

        F_exon_G = nx.Graph()
        R_exon_G = nx.Graph()
        ss_G = nx.DiGraph() 

        all_uniq_exons = {}
        print("n-features:", len(rec.features))
        n_kept = 0
        for feature in rec.features:
            if feature.type in features:
                if len(feature.sub_features)<min_exons:
                    continue
                n_kept +=1
                t = Transcript.init_from_feature(contig, feature)
                #alt_ss = t.get_all_3pSS(contig_seq)

                if t.strand == 1:
                    t.get_splice_junctions(F_3p_5p_ss, 
                                           F_5p_3p_ss, 
                                           F_exon_s_e, 
                                           F_exon_e_s, 
                                           F_5p_to_gene_info, 
                                           F_3p_to_gene_info,
                                           F_le_s_e,
                                           F_le_e_s,
                                           F_fe_s_e,
                                           F_fe_e_s,
                                           all_uniq_exons,
                                           F_LR_intron_interval_tree,
                                           F_added_LR_intron_intervals,
                                           F_exon_G,
                                           ss_G)
                else:
                    t.get_splice_junctions(R_3p_5p_ss, 
                                           R_5p_3p_ss, 
                                           R_exon_s_e, 
                                           R_exon_e_s, 
                                           R_5p_to_gene_info, 
                                           R_3p_to_gene_info,
                                           R_le_s_e,
                                           R_le_e_s,
                                           R_fe_s_e,
                                           R_fe_e_s,
                                           all_uniq_exons,
                                           R_LR_intron_interval_tree,
                                           R_added_LR_intron_intervals,
                                           R_exon_G,
                                           ss_G)
        
        print("n-features kept:", n_kept)
        
        SGs_by_contig[contig] = SpliceGraph(contig, F_3p_5p_ss = F_3p_5p_ss, 
                                                    F_5p_3p_ss = F_5p_3p_ss, 
                                                    R_3p_5p_ss = R_3p_5p_ss, 
                                                    R_5p_3p_ss = R_5p_3p_ss, 
                                                    F_exon_s_e = F_exon_s_e,
                                                    F_exon_e_s = F_exon_e_s,
                                                    R_exon_s_e = R_exon_s_e,
                                                    R_exon_e_s = R_exon_e_s,
                                                    F_5p_to_gene_info = F_5p_to_gene_info,
                                                    F_3p_to_gene_info = F_3p_to_gene_info,
                                                    R_5p_to_gene_info = R_5p_to_gene_info,
                                                    R_3p_to_gene_info = R_3p_to_gene_info,
                                                    F_le_s_e = F_le_s_e,
                                                    F_le_e_s = F_le_e_s,
                                                    R_le_s_e = R_le_s_e,
                                                    R_le_e_s = R_le_e_s,
                                                    F_fe_s_e = F_fe_s_e,
                                                    F_fe_e_s = F_fe_e_s,
                                                    R_fe_s_e = R_fe_s_e,
                                                    R_fe_e_s = R_fe_e_s,
                                                    all_uniq_exons = all_uniq_exons,
                                                    F_LR_intron_interval_tree=F_LR_intron_interval_tree,
                                                    R_LR_intron_interval_tree=R_LR_intron_interval_tree,
                                                    F_exon_G = F_exon_G,
                                                    R_exon_G = R_exon_G,
                                                    ss_G = ss_G)
    return SGs_by_contig
    
#for rec in GFF_parser.parse_in_parts(open(o.fn_gff), limit_info=limit_info):
#    #iterate over chromosomes
#    contig_seq = fa.get_sequence(rec.id)
#    
#    D_to_A, A_to_D = get_donor_acceptor_pairs(rec)
#    for feature in rec.features:
#        if feature.type in ["mRNA", "transcript"]:
#            t = transcript.init_from_feature(feature)
#            alt_ss = t.get_all_3pSS(contig_seq)
#            print alt_ss
