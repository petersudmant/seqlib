import pdb
import itertools
import networkx as nx

from transcript import Transcript

class MisoUtils(object):

    def __init__(self, **kwargs):
        self.sgs_by_contig = kwargs.get("sg_by_contig")
        self.contig_sizes = kwargs.get("contig_sizes")
    
    def make_transcript(self, sg, contig, strand_d, EXONS): 
        EXON_str =  ":".join(["%d,%d"%(e[0],e[1]) for e in EXONS])
        FEATURE_ID = "{contig}:{es}".format(contig=contig, es=EXON_str)
        """ 
        key information off of first exon 5p
        also, note that we are assuming US to DS order of exons
        """
        e1_5p = strand_d=="FWD" and EXONS[0][1] or EXONS[0][0] 
        e1_3p = strand_d=="FWD" and EXONS[0][0] or EXONS[0][1] 
        """
        get gene name based on 5', or if 5' not hashed, because
        you got the SCX, then us the 3'
        """
        if sg.get_gene_name(strand_d, ss_5p=e1_5p):
            GENE_NAME = sg.get_gene_name(strand_d, ss_5p=e1_5p)
            GENE_ID = sg.get_gene_ID(strand_d, ss_5p=e1_5p)
        else:
            GENE_NAME = sg.get_gene_name(strand_d, ss_3p=e1_3p)
            GENE_ID = sg.get_gene_ID(strand_d, ss_3p=e1_3p)

        G_START = min([e[0] for e in EXONS])
        G_END = max([e[1] for e in EXONS])
        STRAND = strand_d == "FWD" and 1 or -1
        
        assert G_START<G_END, "transcript start > transcript end"

        t = Transcript(contig = contig, 
                       feature_ID = FEATURE_ID,
                       exons = EXONS, 
                       gene_name = GENE_NAME,
                       gene_ID = GENE_ID,
                       g_start = G_START,
                       g_end = G_END,
                       strand = STRAND)
        return t 

    def get_SE(self, sg, d_5p, a_3ps, strand):
        if strand=="REV":
            a_3ps = sorted(a_3ps, reverse=strand=="REV")
        else:
            a_3ps = sorted(a_3ps)
        
        i_5p_3p, i_3p_5p, e_5p_3p, e_3p_5p = sg.get_intron_exon_juncs(strand)

        us_exon = sg.get_csx(d_5p, strand, ss_type_5p=True)
                                             
        SEs = []
        n_3ps = len(a_3ps)
        for i in xrange(n_3ps): 
            for j in xrange(i, n_3ps): 
                if i==j: continue
                us_3p = a_3ps[i]
                ds_3p = a_3ps[j]
                us_5ps = set(e_3p_5p[us_3p])
                ds_5ps = set(i_3p_5p[ds_3p])
                shared = us_5ps.intersection(ds_5ps)
                if len(shared)>0:
                    ds_exon = sg.get_csx(ds_3p, strand, ss_type_3p=True)
                    for shared_5p in shared: 
                        SEs.append([us_exon, sorted([us_3p,shared_5p]), ds_exon])
        return SEs
    
    def define_SE_events(self, fn_out_gff, fn_out_bed, source="annot"):
        F_gff = open(fn_out_gff,'w')
        F_bed = open(fn_out_bed,'w')

        for contig, sg in self.sgs_by_contig.items():
            F_R_ss_juncs = { "FWD" : sg.F_5p_3p_ss, 
                             "REV" : sg.R_5p_3p_ss }
            for strand_d, ss_juncs in F_R_ss_juncs.items():
                for d_5p, a_3ps in ss_juncs.iteritems(): 
                    SEs = self.get_SE(sg, d_5p, a_3ps, strand_d)
                    for SE in SEs:
                        us_exon, alt_exon, ds_exon = SE
                        exon_paths = {"A":[0,2], "B":[0,1,2]}
                        EXONS  = [us_exon, alt_exon, ds_exon]
                        trans =  self.make_transcript(sg, contig, strand_d, EXONS) 
                        
                        gff_s = trans.gff_string(exon_paths, source)
                        bed_s = trans.bed_string(exon_paths, source)
                        F_gff.write(gff_s)
                        F_bed.write(bed_s)

    def get_A3SS(self, sg, d_5p, a_3ps, strand):
        """
        test if any of the a_3ps share exon ends 
        """
        if strand=="REV":
            a_3ps = sorted(a_3ps, reverse=strand=="REV")
        else:
            a_3ps = sorted(a_3ps)
        
        i_5p_3p, i_3p_5p, e_5p_3p, e_3p_5p = sg.get_intron_exon_juncs(strand)
        us_exon = sg.get_csx(d_5p, strand, ss_type_5p=True)

        A3SSs = []
        n_3ps = len(a_3ps)
        for i in xrange(n_3ps): 
            for j in xrange(i, n_3ps): 
                if i==j: continue
                us_3p = a_3ps[i]
                ds_3p = a_3ps[j]
                us_5ps = set(e_3p_5p[us_3p])
                ds_5ps = set(e_3p_5p[ds_3p])
                shared = us_5ps.intersection(ds_5ps)
                if len(shared)>0:
                    #don't need to iterate through all shared because 
                    #taking common shortest anyway
                    ds_exon = sg.get_csx(ds_3p, strand, ss_type_3p=True)
                    ds_exon_5p = strand=="FWD" and ds_exon[1] or ds_exon[0]
                    A3SSs.append([us_exon, sorted([us_3p,ds_exon_5p]), ds_exon])
        return A3SSs
    
    def define_A3SS_events(self, fn_out_gff, fn_out_bed, source="annot"):
        F_gff = open(fn_out_gff,'w')
        F_bed = open(fn_out_bed,'w')

        for contig, sg in self.sgs_by_contig.items():
            F_R_ss_juncs = { "FWD" : sg.F_5p_3p_ss, 
                             "REV" : sg.R_5p_3p_ss }
            for strand_d, ss_juncs in F_R_ss_juncs.items():
                for d_5p, a_3ps in ss_juncs.iteritems(): 
                    
                    A3SSs = self.get_A3SS(sg, d_5p, a_3ps, strand_d)
                    for A3SS in A3SSs:
                        us_exon, alt_us_exon, alt_ds_exon = A3SS
                        exon_paths = {"A":[0,1], "B":[0,2]}
                        EXONS  = [us_exon, alt_us_exon, alt_ds_exon]
                        trans =  self.make_transcript(sg, contig, strand_d, EXONS)
                        
                        gff_s = trans.gff_string(exon_paths, source)
                        bed_s = trans.bed_string(exon_paths, source)
                        F_gff.write(gff_s)
                        F_bed.write(bed_s)

    def get_A5SS(self, sg, a_3p, d_5ps, strand):
        """
        test if any of the a_3ps share exon ends 
        """
        if strand=="REV":
            d_5ps = sorted(d_5ps, reverse=strand=="REV")
        else:
            d_5ps = sorted(d_5ps)
        
        i_5p_3p, i_3p_5p, e_5p_3p, e_3p_5p = sg.get_intron_exon_juncs(strand)
        ds_exon = sg.get_csx(a_3p, strand, ss_type_3p=True)

        A5SSs = []
        n_5ps = len(d_5ps)
        for i in xrange(n_5ps): 
            for j in xrange(i, n_5ps): 
                if i==j: continue
                us_5p = d_5ps[i]
                ds_5p = d_5ps[j]
                us_3ps = set(e_5p_3p[us_5p])
                ds_3ps = set(e_5p_3p[ds_5p])
                shared = us_3ps.intersection(ds_3ps)
                if len(shared)>0:
                    #don't need to iterate through all shared because 
                    #taking common shortest anyway
                    us_exon = sg.get_csx(us_5p, strand, ss_type_5p=True)
                    us_exon_3p = strand=="FWD" and us_exon[0] or us_exon[1]
                    A5SSs.append([us_exon, sorted([us_exon_3p,ds_5p]), ds_exon])
        return A5SSs
        
    def define_A5SS_events(self, fn_out_gff, fn_out_bed, source="annot"):
        F_gff = open(fn_out_gff,'w')
        F_bed = open(fn_out_bed,'w')

        for contig, sg in self.sgs_by_contig.items():
            F_R_ss_juncs = { "FWD" : sg.F_3p_5p_ss, 
                             "REV" : sg.R_3p_5p_ss }
            for strand_d, ss_juncs in F_R_ss_juncs.items():
                for a_3p, d_5ps in ss_juncs.iteritems(): 
                    
                    A5SSs = self.get_A5SS(sg, a_3p, d_5ps, strand_d)
                    for A5SS in A5SSs:
                        alt_us_exon, alt_ds_exon, ds_exon = A5SS
                        exon_paths = {"A":[0,2], "B":[1,2]}
                        EXONS  = [alt_us_exon, alt_ds_exon, ds_exon]
                        trans =  self.make_transcript(sg, contig, strand_d, EXONS)
                        
                        gff_s = trans.gff_string(exon_paths, source)
                        bed_s = trans.bed_string(exon_paths, source)
                        F_gff.write(gff_s)
                        F_bed.write(bed_s)
    
    def get_MXE(self, sg, d_5p, a_3ps, strand):
        if strand=="REV":
            a_3ps = sorted(a_3ps, reverse=strand=="REV")
        else:
            a_3ps = sorted(a_3ps)
        
        i_5p_3p, i_3p_5p, e_5p_3p, e_3p_5p = sg.get_intron_exon_juncs(strand)
        le_s_e, le_e_s, fe_s_e, fe_e_s = sg.get_first_last_exons(strand)

        us_exon = sg.get_csx(d_5p, strand, ss_type_5p=True)
                                             
        MXEs = []
        n_3ps = len(a_3ps)
        for i,j in itertools.combinations(range(n_3ps),2):
            assert i!=j
            us_alt_3p = a_3ps[i]
            ds_alt_3p = a_3ps[j]

            us_alt_5ps = set(e_3p_5p[us_alt_3p])
            ds_alt_5ps = set(e_3p_5p[ds_alt_3p])
            shared_alt_5ps = us_alt_5ps.intersection(ds_alt_5ps)
            """
            make sure no exons ending in the same place, 
            otherwise, simply an alt 3' ss
            AND
            make sure the two alts don't splice to each other
            """
            ds_alt_3p_5pdonors = set(i_3p_5p[ds_alt_3p])
            ds_alt_3p_5pdonors_shared = ds_alt_3p_5pdonors.intersection(us_alt_5ps)

            if len(shared_alt_5ps)==0 and len(ds_alt_3p_5pdonors_shared)==0:
                for us_alt_5p in list(us_alt_5ps):
                    for ds_alt_5p in list(ds_alt_5ps):
                        """
                        if last exon, don't try o go further
                        """
                        if us_alt_5p in le_s_e or us_alt_5p in le_e_s: continue
                        if ds_alt_5p in le_s_e or ds_alt_5p in le_e_s: continue

                        us_le_3ps = set(i_5p_3p[us_alt_5p])
                        ds_le_3ps = set(i_5p_3p[ds_alt_5p])
                        
                        shared_le_3ps = us_le_3ps.intersection(ds_le_3ps)
                        for shared_le_3p in shared_le_3ps:
                            ds_exon = sg.get_csx(shared_le_3p, strand, ss_type_3p=True)
                            alt_us_exon = sorted([us_alt_3p,us_alt_5p])
                            alt_ds_exon = sorted([ds_alt_3p,ds_alt_5p])
                            MXEs.append([us_exon, alt_us_exon, alt_ds_exon, ds_exon])
        return MXEs
    
    def define_MXE_events(self, fn_out_gff, fn_out_bed, source="annot"):
        F_gff = open(fn_out_gff,'w')
        F_bed = open(fn_out_bed,'w')

        for contig, sg in self.sgs_by_contig.items():
            F_R_ss_juncs = { "FWD" : sg.F_5p_3p_ss, 
                             "REV" : sg.R_5p_3p_ss }
            for strand_d, ss_juncs in F_R_ss_juncs.items():
                for d_5p, a_3ps in ss_juncs.iteritems(): 
                    MXEs = self.get_MXE(sg, d_5p, a_3ps, strand_d)
                    for MXE in MXEs:
                        us_exon, alt_us_exon, alt_ds_exon, ds_exon = MXE
                        exon_paths = {"A":[0,1,3], "B":[0,2,3]}
                        EXONS  = [us_exon, alt_us_exon, alt_ds_exon, ds_exon]
                        trans =  self.make_transcript(sg, contig, strand_d, EXONS)
                        
                        gff_s = trans.gff_string(exon_paths, source)
                        bed_s = trans.bed_string(exon_paths, source)
                        F_gff.write(gff_s)
                        F_bed.write(bed_s)

    def overlap(self, ex1, ex2):
        #exons in s,e format 
        return (ex1[0]<=ex2[1] and ex2[0]<=ex1[1])
    
    def get_AFE(self, sg, connected_exs, strand):
        le_s_e, le_e_s, fe_s_e, fe_e_s = sg.get_first_last_exons(strand)
        i_5p_3p, i_3p_5p, e_5p_3p, e_3p_5p = sg.get_intron_exon_juncs(strand)

        curr_fexs = [ex for ex in connected_exs if ex[0] in fe_s_e]
        
        if strand=="REV":
            curr_fexs = sorted(curr_fexs, reverse=True)
            ex_to_donor = lambda x: x[0]
        else:
            curr_fexs = sorted(curr_fexs)
            ex_to_donor = lambda x: x[1]
        
        ##now collapse fexs to non-overlapping
        fex_uniq_5ps = list(set([ex_to_donor(fex) for fex in curr_fexs]))
        csx_fexs = [sg.get_csx(d_5p, strand, ss_type_5p=True) for d_5p in fex_uniq_5ps]
        AFEs = []
        n_fexs = len(csx_fexs)
        
        for i,j in itertools.combinations(range(n_fexs),2):
            assert i!=j
            us_fex = csx_fexs[i] #distal, us
            ds_fex = csx_fexs[j] #proximal, ds
            if self.overlap(us_fex,ds_fex): continue
            AFEs.append([us_fex,ds_fex])
        return AFEs

    def __define_AFE_events(self, fn_out_gff, fn_out_bed, source="annot"):
        F_gff = open(fn_out_gff,'w')
        F_bed = open(fn_out_bed,'w')

        for contig, sg in self.sgs_by_contig.items():
            F_R_ex_graphs = { "FWD" : sg.F_exon_G, 
                              "REV" : sg.R_exon_G }
            for strand_d, ex_G in F_R_ex_graphs.items():
                for connected_exs in nx.connected_components(ex_G): 
                    AFEs = self.get_AFE(sg, connected_exs, strand_d)
                    for AFE in AFEs:
                        EXONS = AFE
                        exon_paths = {"A":[0], "B":[1]}
                        trans =  self.make_transcript(sg, contig, strand_d, EXONS)
                        gff_s = trans.gff_string(exon_paths, source)
                        bed_s = trans.bed_string(exon_paths, source)
                        F_gff.write(gff_s)
                        F_bed.write(bed_s)

    def get_ALE(self, sg, connected_exs, strand):
        """
        TO DO _ remove overlapping exons..., just take the... shortest?
        example: is this really an ALE?
        chr19:46,325,214-46,335,433
        """
        le_s_e, le_e_s, fe_s_e, fe_e_s = sg.get_first_last_exons(strand)
        i_5p_3p, i_3p_5p, e_5p_3p, e_3p_5p = sg.get_intron_exon_juncs(strand)

        curr_lexs = [ex for ex in connected_exs if ex[0] in le_s_e]
        
        if strand=="REV":
            curr_lexs = sorted(curr_lexs, reverse=True)
            ex_to_acceptor = lambda x: x[1]
        else:
            curr_lexs = sorted(curr_lexs)
            ex_to_acceptor = lambda x: x[0]
        
        ##now collapse lexs to non-overlapping
        lex_uniq_3ps = list(set([ex_to_acceptor(lex) for lex in curr_lexs]))
        csx_lexs = [sg.get_csx(a_3p, strand, ss_type_3p=True) for a_3p in lex_uniq_3ps]
        ALEs = []
        n_lexs = len(csx_lexs)
        
        for i,j in itertools.combinations(range(n_lexs),2):
            assert i!=j
            us_lex = csx_lexs[i] #distal, us
            ds_lex = csx_lexs[j] #proximal, ds
            if self.overlap(us_lex,ds_lex): continue
            ALEs.append([us_lex,ds_lex])
        return ALEs
    
        
    def __define_ALE_events(self, fn_out_gff, fn_out_bed, source="annot"):
        F_gff = open(fn_out_gff,'w')
        F_bed = open(fn_out_bed,'w')

        for contig, sg in self.sgs_by_contig.items():
            F_R_ex_graphs = { "FWD" : sg.F_exon_G, 
                              "REV" : sg.R_exon_G }
            for strand_d, ex_G in F_R_ex_graphs.items():
                for connected_exs in nx.connected_components(ex_G): 
                    ALEs = self.get_ALE(sg, connected_exs, strand_d)
                    for ALE in ALEs:
                        EXONS = ALE
                        exon_paths = {"A":[0], "B":[1]}
                        trans =  self.make_transcript(sg, contig, strand_d, EXONS)
                        gff_s = trans.gff_string(exon_paths, source)
                        bed_s = trans.bed_string(exon_paths, source)
                        F_gff.write(gff_s)
                        F_bed.write(bed_s)


    def define_AFE_events(self, fn_out_gff, fn_out_bed, source="annot"):
        self.define_simple_AFE_events(fn_out_gff, fn_out_bed, source=source)

    def define_ALE_events(self, fn_out_gff, fn_out_bed, source="annot"):
        self.define_simple_ALE_events(fn_out_gff, fn_out_bed, source=source)
        
    """collapse a list of exons into their unique components"""
    def collapse_exons(self, exs):
        exs = sorted(exs)
        collapsed = [exs[0]]
        
        for ex in exs:
            assert ex[0]>=collapsed[-1][0]
            
            if ex[0]<collapsed[-1][1]:
                collapsed[-1] = tuple([collapsed[-1][0], ex[1]])
            else:
                collapsed.append(ex)
        
        return collapsed

    def get_simple_ALE(self, sg, connected_exs, strand):
        le_s_e, le_e_s, fe_s_e, fe_e_s = sg.get_first_last_exons(strand)
        i_5p_3p, i_3p_5p, e_5p_3p, e_3p_5p = sg.get_intron_exon_juncs(strand)

        curr_lexs = [ex for ex in connected_exs if ex[0] in le_s_e]
        
        if strand=="REV":
            curr_lexs = sorted(curr_lexs, reverse=True)
            ex_to_acceptor = lambda x: x[1]
        else:
            curr_lexs = sorted(curr_lexs)
            ex_to_acceptor = lambda x: x[0]
        
        lex_uniq_3ps = list(set([ex_to_acceptor(lex) for lex in curr_lexs]))
        csx_lexs = [sg.get_csx(a_3p, strand, ss_type_3p=True) for a_3p in lex_uniq_3ps]

        ##now collapse lexs to non-overlapping
        csx_lexs = self.collapse_exons(csx_lexs)
        
        """could exclude anything with only 1"""
        ALEs = csx_lexs
        
        return ALEs

    """
    simple ale/afes are only the first/last exon in question
    as done my jason... 
    """
    def define_simple_ALE_events(self, fn_out_gff, fn_out_bed, source="annot"):
        F_gff = open(fn_out_gff,'w')
        F_bed = open(fn_out_bed,'w')

        for contig, sg in self.sgs_by_contig.items():
            F_R_ex_graphs = { "FWD" : sg.F_exon_G, 
                              "REV" : sg.R_exon_G }
            for strand_d, ex_G in F_R_ex_graphs.items():
                for connected_exs in nx.connected_components(ex_G): 
                    ALEs = self.get_simple_ALE(sg, connected_exs, strand_d)
                    EXONS = ALEs
                    exon_paths = {chr(i + ord('A')):[i] for i,ex in enumerate(EXONS)}
                    trans =  self.make_transcript(sg, contig, strand_d, EXONS)
                    gff_s = trans.gff_string(exon_paths, source)
                    bed_s = trans.bed_string(exon_paths, source)
                    F_gff.write(gff_s)
                    F_bed.write(bed_s)

    def get_simple_AFE(self, sg, connected_exs, strand):
        le_s_e, le_e_s, fe_s_e, fe_e_s = sg.get_first_last_exons(strand)
        i_5p_3p, i_3p_5p, e_5p_3p, e_3p_5p = sg.get_intron_exon_juncs(strand)

        curr_fexs = [ex for ex in connected_exs if ex[0] in fe_s_e]
        
        if strand=="REV":
            curr_fexs = sorted(curr_fexs, reverse=True)
            ex_to_donor = lambda x: x[0]
        else:
            curr_fexs = sorted(curr_fexs)
            ex_to_donor = lambda x: x[1]
        
        ##now collapse fexs to non-overlapping
        fex_uniq_5ps = list(set([ex_to_donor(fex) for fex in curr_fexs]))
        csx_fexs = [sg.get_csx(d_5p, strand, ss_type_5p=True) for d_5p in fex_uniq_5ps]

        ##now collapse fexs to non-overlapping
        csx_fexs = self.collapse_exons(csx_fexs)
        
        """could exclude anything with only 1"""
        AFEs = csx_fexs
        
        return AFEs
        
    def define_simple_AFE_events(self, fn_out_gff, fn_out_bed, source="annot"):
        F_gff = open(fn_out_gff,'w')
        F_bed = open(fn_out_bed,'w')

        for contig, sg in self.sgs_by_contig.items():
            F_R_ex_graphs = { "FWD" : sg.F_exon_G, 
                              "REV" : sg.R_exon_G }
            for strand_d, ex_G in F_R_ex_graphs.items():
                for connected_exs in nx.connected_components(ex_G): 
                    AFEs = self.get_simple_AFE(sg, connected_exs, strand_d)
                    EXONS = AFEs
                    exon_paths = {chr(i + ord('A')):[i] for i,ex in enumerate(EXONS)}
                    trans =  self.make_transcript(sg, contig, strand_d, EXONS)
                    gff_s = trans.gff_string(exon_paths, source)
                    bed_s = trans.bed_string(exon_paths, source)
                    F_gff.write(gff_s)
                    F_bed.write(bed_s)

    def get_RIs(self, sg, a_3p, d_5ps, strand):
        if strand=="REV":
            d_5ps = sorted(d_5ps, reverse=strand=="REV")
        else:
            d_5ps = sorted(d_5ps)
        
        i_5p_3p, i_3p_5p, e_5p_3p, e_3p_5p = sg.get_intron_exon_juncs(strand)
        
        ds_ex_5ps = set(e_3p_5p[a_3p])

        RIs = []
        for us_5p in d_5ps:
            for us_ex_3p in e_5p_3p[us_5p]:
                us_ex_5ps = set(e_3p_5p[us_ex_3p])
                intersection = us_ex_5ps.intersection(ds_ex_5ps)
                for shared_5p in intersection:
                    RI_exon = tuple(sorted([us_ex_3p, shared_5p]))
                    us_ex = tuple(sorted([us_ex_3p,us_5p]))
                    ds_ex = tuple(sorted([a_3p,shared_5p]))
                    RIs.append([us_ex, ds_ex, RI_exon])
        return RIs


    def define_RI_events(self, fn_out_gff, fn_out_bed, source="annot"):
        F_gff = open(fn_out_gff,'w')
        F_bed = open(fn_out_bed,'w')

        for contig, sg in self.sgs_by_contig.items():
            F_R_ss_juncs = { "FWD" : sg.F_3p_5p_ss, 
                             "REV" : sg.R_3p_5p_ss }
            for strand_d, ss_juncs in F_R_ss_juncs.items():
                for a_3p, d_5ps in ss_juncs.iteritems(): 
                    RIs = self.get_RIs(sg, a_3p, d_5ps, strand_d)
                    for RI in RIs:
                        EXONS = RI #us_exon, ds_ex, RI_exon
                        exon_paths = {"A":[0,1], "B":[2]}
                        trans =  self.make_transcript(sg, contig, strand_d, EXONS)
                        gff_s = trans.gff_string(exon_paths, source)
                        bed_s = trans.bed_string(exon_paths, source)
                        F_gff.write(gff_s)
                        F_bed.write(bed_s)



