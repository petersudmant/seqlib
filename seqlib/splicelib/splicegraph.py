import numpy as np
import pdb
import argparse
import BCBio.GFF as GFF
import re
import string 

from transcript import Transcript
from collections import defaultdict
from sys import stderr
from fastahack import FastaHack

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
        
        self.F_3p_5p_ss = kwargs["F_3p_5p_ss"]
        self.F_5p_3p_ss = kwargs["F_5p_3p_ss"]
        self.R_3p_5p_ss = kwargs["R_3p_5p_ss"]
        self.R_5p_3p_ss = kwargs["R_5p_3p_ss"]
        
        self.F_exon_s_e = kwargs["F_exon_s_e"]
        self.F_exon_e_s = kwargs["F_exon_e_s"]
        self.R_exon_s_e = kwargs["R_exon_s_e"]
        self.R_exon_e_s = kwargs["R_exon_e_s"]
        
        self.F_5p_to_gene_info = kwargs["F_5p_to_gene_info"]
        self.F_3p_to_gene_info = kwargs["F_3p_to_gene_info"]
        self.R_5p_to_gene_info = kwargs["R_5p_to_gene_info"]
        self.R_3p_to_gene_info = kwargs["R_3p_to_gene_info"]
        
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
    
    

    def get_micro_exon(self, seq, F_gff, F_bed, n=33):
        """
        define micro_exon

        SD1, SA1, SD2, SA2 where 
        SA1 and SA2 are within n bp of each other

        again, assuming common shortest exon here
        """

        """
        get all 5' 3' ss pairs and search in-between for lil' guys
        """
        for strand_d, _5p_3p_ss in { "FWD" : self.F_5p_3p_ss, "REV" : self.R_5p_3p_ss }.iteritems():
            
            if strand_d == "FWD":
                _5p_to_gene_inf = self.F_5p_to_gene_info
                _3p_to_gene_inf = self.F_3p_to_gene_info
            else:
                _5p_to_gene_inf = self.R_5p_to_gene_info
                _3p_to_gene_inf = self.R_3p_to_gene_info

            for ss_5p, ss_3p_list in _5p_3p_ss.iteritems():
                
                us_exon = self.get_common_shortest_exon(ss_5p, strand_d, ss_type_5p=True)
                us_gene_inf = _5p_to_gene_info[ss_5p]
                
                for ss_3p in ss_3p_list:
                    ds_exon = self.get_common_shortest_exon(ss_3p,  strand_d, ss_type_3p=True)
                    ds_gene_inf = _3p_to_gene_info[ss_3p]
                    
                    ss_5p_seq = strand_d == "FWD" and "GT" or "AC"
                    ss_3p_seq = strand_d == "FWD" and "AG" or "CT"
                    
                    s,e = sorted([ss_5p, ss_3p])
                    s,e = s+2,e-2

                    alt_5p_ss = np.array([m.start()+s for m in re.finditer(ss_5p_seq, seq[s:e].upper())])
                    alt_3p_ss = np.array([m.start()+s for m in re.finditer(ss_3p_seq, seq[s:e].upper())])

                    #OK, now ID the chunks that fit the bill!, ie, ss- < n bases -ss



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
                for alt_annot_5_or_3 in alt_annot_ss_5_or_3_list: 
                    if get_3p:
                        curr_3p = alt_annot_5_or_3
                        curr_5p = fixed_ss_3_or_5
                    if get_5p:
                        curr_3p = fixed_ss_3_or_5   
                        curr_5p = alt_annot_5_or_3

                ds_exon = self.get_common_shortest_exon(curr_3p, strand_d, ss_type_3p=True)
                us_exon = self.get_common_shortest_exon(curr_5p, strand_d, ss_type_5p=True)
                
                """
                GET GENE INFO
                """
                if strand_d == "FWD":
                    us_gene_inf = self.F_5p_to_gene_info[curr_5p]
                    ds_gene_inf = self.F_3p_to_gene_info[curr_3p]
                else:
                    us_gene_inf = self.R_5p_to_gene_info[curr_5p]
                    ds_gene_inf = self.R_3p_to_gene_info[curr_3p]

                gene_inf = []
                for inf in us_gene_inf+ds_gene_inf:
                    if not inf in gene_inf: gene_inf.append(inf)
            
                """
                k is the number of bp up/down of the annotated 
                5p/3p ss to search for now ones, ie, don't search past 
                annotated OTHER 5' splice site!
                """
                
                if get_3p:
                    ss_seq = strand_d == "FWD" and "AG" or "CT"
                else:
                    ss_seq = strand_d == "FWD" and "GT" or "AC"
                        
                max_k = min(n, abs(curr_5p-curr_3p))
            
                min_alt, max_alt = alt_annot_5_or_3-k, alt_annot_5_or_3+k
                if strand_d == "FWD":
                    #seq_s, seq_e = annot_3p_ss-k, annot_3p_ss-2
                    if get_5p:
                        seq_s, seq_e = max(min_alt, us_exon[0]+1), min(max_alt, ds_exon[0])
                    else: 
                        seq_s, seq_e = max(min_alt, us_exon[1]+1), min(max_alt, ds_exon[1])
                else:
                    if get_5p:
                        seq_s, seq_e = max(min_alt, ds_exon[1]+1), min(max_alt, us_exon[1])
                    else:
                        seq_s, seq_e = max(min_alt, ss_exon[0]+1), min(max_alt, us_exon[0])
                """
                offset by a 2 if in fwd dir as you want exon 
                EXON start, not the, position of the ss
                for rev, the position is fine already
                """
                delta = strand_d == "FWD" and 2 or 0 
                alt_ss_list = np.array([m.start()+seq_s+delta for m in re.finditer(ss_seq, seq[seq_s:seq_e].upper())])
                #############
                for alt_ss in alt_ss_list:
                    if alt_ss == alt_annot_5_or_3: 
                        continue

                    if alt_ss in alt_annot_ss_5_or_3_list: 
                        source = "annotated"
                    else:
                        source = "novel"
                    
                    exon_paths = {"A":[0,1], "B":[0,2]}
                    if "FWD" and get_3p:
                        alt_exon=[alt_ss, ds_exon[1]]
                    elif "REV" and get_3p:
                        alt_exon=[ds_exon[0], alt_ss]
                    elif "FWD" and get_5p:
                        alt_exon=[us_exon[0], alt_ss]
                    elif "REV" and get_5p:
                        alt_exon=[alt_ss, us_exon[1]]

                    EXONS  = [us_exon, alt_exon, ds_exon]
                    FEATURE_ID = "_".join(["%s:%d-%d"%(self.contig, e[0], e[1]) for e in EXONS])
                    GENE_NAME = ",".join(["%s"%gi['gene_name'] for gi in gene_inf])
                    GENE_ID = "_".join(["%s"%gi['gene_ID'] for gi in gene_inf])
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

                    gff_s = alt_T.gff_string(exon_paths, source)
                    gff_novel_s = alt_T.gff_string({"A":[0,1]}, source)

                    bed_s = alt_T.bed_string(exon_paths, source, True)
                    junc_tups = alt_T.junc_tuples(exon_paths, source, True)
                    F_gff.write(gff_s)
                    F_novel_gff.write(gff_novel_s)
                    F_bed.write(bed_s)
                    junc_writer.write(junc_tups)


    def __DEPRECATED__get_NAGNAGs(self, seq, F_gff, F_novel_gff, F_bed, junc_writer, n=60):
        """
        define NAGNAG
        SD SA1,SA2 where SA1 and SA2 are within n bp of each other
        """
        
        for strand_d, _5p_3p_ss in { "FWD" : self.F_5p_3p_ss, "REV" : self.R_5p_3p_ss }.iteritems():
            for ss_5p, ss_3p_list in _5p_3p_ss.iteritems():
                annotated_ss_3p = np.sort(np.array(ss_3p_list))

                max_ss_3p = annotated_ss_3p[-1]
                min_ss_3p = annotated_ss_3p[0]
                
                if strand_d == "FWD":
                    annot_3p_ss = max_ss_3p
                else:
                    annot_3p_ss = min_ss_3p

                ds_exon = self.get_common_shortest_exon(annot_3p_ss, strand_d, ss_type_3p=True)
                us_exon = self.get_common_shortest_exon(ss_5p,  strand_d, ss_type_5p=True)
                
                """
                GET GENE INFO
                """
                if strand_d == "FWD":
                    us_gene_inf = self.F_5p_to_gene_info[ss_5p]
                    ds_gene_inf = self.F_3p_to_gene_info[annot_3p_ss]
                else:
                    us_gene_inf = self.R_5p_to_gene_info[ss_5p]
                    ds_gene_inf = self.R_3p_to_gene_info[annot_3p_ss]

                gene_inf = []
                for inf in us_gene_inf+ds_gene_inf:
                    if not inf in gene_inf: gene_inf.append(inf)
            
                """
                k is the number of bp up/down of the annotated 
                3p ss to search for now ones, ie, don't search past 
                annotated 5' splice site!
                """
                
                k = min(n, abs(annot_3p_ss-ss_5p))
                ss_seq = strand_d == "FWD" and "AG" or "CT"
                if strand_d == "FWD":
                    #seq_s, seq_e = annot_3p_ss-k, annot_3p_ss-2
                    seq_s, seq_e = annot_3p_ss-k, min(annot_3p_ss+k+2, ds_exon[1]-1)
                else:
                    #seq_s, seq_e = annot_3p_ss+2, annot_3p_ss+k
                    seq_s, seq_e = max(ds_exon[0]+1, annot_3p_ss-k), annot_3p_ss+k
                """
                offset by a 2 if in fwd dir as you want exon 
                EXON start, not the, position of the ss
                for rev, the position is fine already
                """
                delta = strand_d == "FWD" and 2 or 0 
                alt_3p_ss = np.array([m.start()+seq_s+delta for m in re.finditer(ss_seq, seq[seq_s:seq_e].upper())])

                for ss_3p in alt_3p_ss:
                    if ss_3p == annot_3p_ss: 
                        continue

                    if ss_3p in annotated_ss_3p: 
                        source = "NAGNNAG_annotated"
                    else:
                        source = "NAGNNAG"
                    
                    exon_paths = {"A":[0,1], "B":[0,2]}
                    NAG_exon = strand_d == "FWD" and [ss_3p, ds_exon[1]] or [ds_exon[0], ss_3p]
                    EXONS  = [us_exon, NAG_exon, ds_exon]
                    FEATURE_ID = "_".join(["%s:%d-%d"%(self.contig, e[0], e[1]) for e in EXONS])
                    GENE_NAME = ",".join(["%s"%gi['gene_name'] for gi in gene_inf])
                    GENE_ID = "_".join(["%s"%gi['gene_ID'] for gi in gene_inf])
                    G_START = min(us_exon[0], ds_exon[0])
                    G_END = max(us_exon[1], ds_exon[1])
                    STRAND = strand_d == "FWD" and 1 or -1
                    
                    nagNnag_T = Transcript(contig = self.contig, 
                                           feature_ID = FEATURE_ID,
                                           exons = EXONS, 
                                           gene_name = GENE_NAME,
                                           gene_ID = GENE_ID,
                                           g_start = G_START,
                                           g_end = G_END,
                                           strand = STRAND)

                    gff_s = nagNnag_T.gff_string(exon_paths, source)
                    gff_novel_s = nagNnag_T.gff_string({"A":[0,1]}, source)

                    bed_s = nagNnag_T.bed_string(exon_paths, source, True)
                    junc_tups = nagNnag_T.junc_tuples(exon_paths, source, True)
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

        F_5p_to_gene_info = {}
        F_3p_to_gene_info = {}
        R_5p_to_gene_info = {}
        R_3p_to_gene_info = {}

        for feature in rec.features:
            if feature.type in features:
                t = Transcript.init_from_feature(contig, feature)
                #alt_ss = t.get_all_3pSS(contig_seq)
                if t.strand == 1:
                    t.get_splice_junctions(F_3p_5p_ss, 
                                           F_5p_3p_ss, 
                                           F_exon_s_e, 
                                           F_exon_e_s, 
                                           F_5p_to_gene_info, 
                                           F_3p_to_gene_info)
                else:
                    t.get_splice_junctions(R_3p_5p_ss, 
                                           R_5p_3p_ss, 
                                           R_exon_s_e, 
                                           R_exon_e_s, 
                                           R_5p_to_gene_info, 
                                           R_3p_to_gene_info)
        
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
                                                    R_3p_to_gene_info = R_3p_to_gene_info)
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
