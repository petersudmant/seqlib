import pdb
from transcript import Transcript

class MisoUtils(object):

    def __init__(self, **kwargs):
        self.sgs_by_contig = kwargs.get("sg_by_contig")
        self.contig_sizes = kwargs.get("contig_sizes")
    
    def make_transcript(self, sg, contig, strand_d, exons, path, d_5p): 
        mx = self.contig_sizes[contig]
        EXONS = [[e[0],min(e[1],mx)] for e in exons]
        
        EXON_str =  ":".join(["%d,%d"%(e[0],e[1]) for e in EXONS])
        FEATURE_ID = "{contig}:{es}".format(contig=contig, es=EXON_str)

        GENE_NAME = sg.get_gene_name(strand_d, ss_5p=d_5p)
        GENE_ID = sg.get_gene_ID(strand_d, ss_5p=d_5p)

        G_START = min([e[0] for e in EXONS])
        G_END = max([e[1] for e in EXONS])
        STRAND = strand_d == "FWD" and 1 or -1
    
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
                        SEs.append([us_exon, [us_3p,shared_5p], ds_exon])
        return SEs
    
    def define_SE_events(self, fn_out, source="annot"):
        F_gff = open(fn_out,'w')
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

                        trans =  self.make_transcript(sg, contig, strand_d, EXONS, exon_paths, d_5p) 
                        gff_s = trans.gff_string(exon_paths, source)
                        F_gff.write(gff_s)

    def get_A3SS(self, sg, d_5p, a_3ps, strand_d):
        pass

    def define_A3SS_events(self, fn_out):
        pass
        F_gff = open(fn_out,'w')
        for contig, sg in self.sgs_by_contig.items():
            F_R_ss_juncs = { "FWD" : sg.F_5p_3p_ss, 
                             "REV" : sg.R_5p_3p_ss }
            for strand_d, ss_juncs in F_R_ss_juncs.items():
                for d_5p, a_3ps in ss_juncs.iteritems(): 
                    
                    A3SSs = self.get_A3SS(sg, d_5p, a_3ps, strand_d)
                    for A3SS in A3SS:
                        us_exon, alt_exon, ds_exon = SE
                        exon_paths = {"A":[0,1], "B":[0,2]}
                        EXONS  = [us_exon, alt_exon, ds_exon]
                        
                        trans =  self.make_transcript(sg, contig, strand_d, EXONS, exon_paths, d_5p) 
                        
                        gff_s = trans.gff_string(exon_paths, source)
                        F_gff.write(gff_s)
    


    def define_A5SS_events(self, fn_out):
        pass
    def define_MXE_events(self, fn_out):
        pass
    def define_RI_events(self, fn_out):
        pass
    def define_ALE_events(self, fn_out):
        pass
    def define_AFE_events(self, fn_out):
        pass



