import pdb
from bx.intervals.intersection import Interval, IntervalTree

class Transcript:

    def __init__(self, **kwargs):
        self.feature_ID = kwargs['feature_ID'] 
        self.gene_name = kwargs['gene_name'] 
        self.gene_ID = kwargs['gene_ID'] 
        self.contig = kwargs['contig'] 
        self.g_start = kwargs['g_start'] 
        self.g_end = kwargs['g_end'] 
        self.exons = kwargs['exons']
        self.strand = kwargs['strand']
    
    def get_tuple_hash(self):
        """
        a tuple describing a transcript uniquely to see if it's already been
        output
        """
        exon_starts = [e[0] for e in self.exons]
        exon_ends = [e[1] for e in self.exons]
        return tuple(sorted(exon_starts+exon_ends))

    def print_out(self):
        print self.feature_ID
        print "\t", self.gene_name #important
        print "\t", self.gene_ID  #important
        print "\t", self.g_start
        print "\t", self.g_end
        print "\t", self.exons
        print "\t", self.strand
    
    def get_gene_info(self):
        return {"gene_name":self.gene_name, "gene_ID": self.gene_ID}
    
    def junc_string(self, exon_paths, source, UCSC=False):
        junc_lines = []
        pattern = "{contig}\t{left}\t{right}\t{strand}"
        
        output_tups = []
        for tup in self.junc_tuples(exon_paths, source, UCSC):
            if not tup in output_tups:
                output_tups.append(tup)
                junc_lines.append(pattern.format(contig=tup[0],
                                                 left = tup[1],
                                                 right= tup[2],
                                                 strand = tup[3]))
        return "%s\n"%"\n".join(junc_lines)
    
    def junc_tuples(self, exon_paths, source, UCSC=False):

        junc_tups = []

        strand = self.strand == 1 and "+" or "-"

        formatted_contig = self.contig
        if UCSC and "chr" not in self.contig:
            formatted_contig  = "chr%s"%self.contig

        for curr_ID, e_path in exon_paths.iteritems():
            curr_exons = [self.exons[i] for i in e_path]
            if strand == "-": curr_exons = curr_exons[::-1] 
            """
            NOTE THE -1 here because your exon coordinates are half open, 
            so, exon END should be -1
            contig, left, right, strand
            """
            for i in xrange(len(curr_exons)-1):
                junc_tups.append(tuple([formatted_contig,
                                        curr_exons[i][1]-1,
                                        curr_exons[i+1][0],
                                        strand]))
        return junc_tups



    def bed_string(self, exon_paths, source, UCSC=False, annot_RGB="000,204,102", novel_RGB="64,64,64"):
        bed_lines = []
        RGB = source=="annotated" and annot_RGB or novel_RGB
        pattern = ("{contig}\t{start}\t{end}\t{name}\t"
                   "0\t{strand}\t{start}\t{end}\t{RGB}\t"
                   "{n_exons}\t{exon_sizes}\t{exon_starts}")
                   
        strand = self.strand == 1 and "+" or "-"
        
        formatted_contig = self.contig
        if UCSC and "chr" not in self.contig:
            formatted_contig  = "chr%s"%self.contig
        for curr_ID, e_path in exon_paths.iteritems():
            curr_exons = [self.exons[i] for i in e_path]

            t_start = min(curr_exons[0][0], curr_exons[-1][0])
            t_end = max(curr_exons[0][1], curr_exons[-1][1])

            if strand == "-": curr_exons = curr_exons[::-1] 
            min_p = min([min(e[0],e[1]) for e in curr_exons])
            bed_lines.append(pattern.format(contig=formatted_contig,
                                            start = t_start,
                                            end = t_end,
                                            name="%s,%s"%(self.feature_ID, curr_ID),
                                            strand = strand,
                                            RGB = RGB,
                                            n_exons = len(curr_exons),
                                            exon_sizes = ",".join(["%d"%(e[1]-e[0]) for e in curr_exons]),
                                            exon_starts = ",".join(["%d"%(e[0]-min_p) for e in curr_exons]) ))

        return "%s\n"%("\n".join(bed_lines))

    def gff_string(self, exon_paths, source):
        """
        gene / mRNA / exon / exon
        ****NOTE****
        GFF SPEC STIPULATES 1-based COORDS!
        http://useast.ensembl.org/info/website/upload/gff.html
        ***********

        SO, +1 to the start only, because standard is 0 based HALF open, 
        this is one 1 based fully closed
        """
        gff_lines = []
        pattern = "{contig}\t{source}\t{_type}\t{start}\t{end}\t.\t{strand}\t.\tID={ID};{other}"
        strand = self.strand == 1 and "+" or "-"
        gff_lines.append( pattern.format(contig=self.contig, 
                                         source = source, 
                                         _type = "gene", 
                                         start = self.g_start+1,
                                         end = self.g_end, 
                                         strand = strand,
                                         ID = self.feature_ID,
                                         other = "Name=%s;gene_ID=%s"%(self.gene_name, self.gene_ID)))
        """
        now, iterate over exon combinations
        """
        for curr_ID, e_path in exon_paths.iteritems():
            mRNA_ID = "%s_%s"%(self.feature_ID, curr_ID)
            
            curr_exons = [self.exons[i] for i in e_path]
            t_start = min(curr_exons[0][0], curr_exons[-1][0])
            t_end = max(curr_exons[0][1], curr_exons[-1][1])

            gff_lines.append( pattern.format(contig=self.contig, 
                                             source = source, 
                                             _type = "mRNA", 
                                             start = t_start+1,
                                             end = t_end,
                                             strand = strand,
                                             ID = mRNA_ID, 
                                             other = "Parent=%s"%(self.feature_ID)))
            for exon_i in e_path:
                e_start, e_end = self.exons[exon_i]
                exon_ID = "%s_%d"%(mRNA_ID,exon_i) 
                gff_lines.append( pattern.format(contig=self.contig, 
                                                 source = source, 
                                                 _type = "exon", 
                                                 start = e_start+1,
                                                 end = e_end, 
                                                 strand = strand,
                                                 ID = exon_ID, 
                                                 other = "Parent=%s"%(mRNA_ID)))
        
        return "%s\n"%("\n".join(gff_lines))

    def get_all_3pSS(self, seq):

        for i in xrange(len(self.exons)-1):
            e1_e = self.exons[i][1]
            e2_s = self.exons[i+1][0]
            
            if self.strand == 1:
                alts = np.array([m.start() for m in re.finditer("AG",seq[e1_e:e2_s])])
            else:
                alts = np.array([m.start() for m in re.finditer("CT",seq[e1_e:e2_s])])
            return alts 
    
    def get_splice_junctions(self, _3p_to_5p, 
                                   _5p_to_3p, 
                                   _exon_s_to_e, 
                                   _exon_e_to_s, 
                                   _5p_to_gene_info, 
                                   _3p_to_gene_info, 
                                   _le_s_e,
                                   _le_e_s,
                                   _fe_s_e,
                                   _fe_e_s,
                                   all_uniq_exons,
                                   LR_intron_interval_tree,
                                   added_LR_intron_intervals,
                                   exon_G, 
                                   ss_G):

        """
        side effect - modifies passed in dicts
        5p_to_3p and 3p_to_5p are strand specific defaultdicts
        all exons "starts" are strictly < exon "ends"
        not even for exons all the 5'/3' desigations are intron ss based
                          ___
         ---->         5'/   \ 3' 
         |         -----/     \-----
                   ___
                3'/   \ 5'      <----
            -----/     \-----        |
        """
       

        """
        exon order is dependent on strand
        """
        if self.strand == 1:
            fe_s, fe_e = self.exons[0] 
            le_s, le_e = self.exons[-1] 
        else:
            fe_s, fe_e = self.exons[-1] 
            le_s, le_e = self.exons[0] 

        _fe_s_e[fe_s] = fe_e
        _fe_e_s[fe_e] = fe_s
        _le_s_e[le_s] = le_e
        _le_e_s[le_e] = le_s
        
        gene_inf = self.get_gene_info()
        
        for i in xrange(len(self.exons)-1):

            e1_s, e1_e = self.exons[i]
            e2_s, e2_e = self.exons[i+1]
            """
            assertions here - all exons processed from left to right
            ie e1 < e2 and ei < ei+1
            """
            assert e1_s < e1_e and e2_s < e2_e
            assert e1_e < e2_s #previously whack due to CDS
            
            #add the exons
            e1_tup, e2_tup = tuple([e1_s, e1_e]), tuple([e2_s, e2_e])
            if not e1_s in _exon_s_to_e: _exon_s_to_e[e1_s] = []
            if not e1_e in _exon_e_to_s: _exon_e_to_s[e1_e] = []
            if not e1_tup in all_uniq_exons: all_uniq_exons[e1_tup] = 1
            if not e2_tup in all_uniq_exons: all_uniq_exons[e2_tup] = 1
            
            if not e1_e in _exon_s_to_e[e1_s]:
                _exon_s_to_e[e1_s].append(e1_e)
            if not e1_s in _exon_e_to_s[e1_e]:
                _exon_e_to_s[e1_e].append(e1_s)
            
            if not tuple([e1_e, e2_s]) in added_LR_intron_intervals:
                added_LR_intron_intervals[tuple([e1_e, e2_s])] = 1
                LR_intron_interval_tree.insert_interval(Interval(e1_e, e2_s))
            
            #add exons to exon graph
            exon_G.add_edge(e1_tup, e2_tup)

            #add the introns
            if self.strand == 1:
                #FWD
                if not e1_e in _5p_to_3p: _5p_to_3p[e1_e] = []
                if not e2_s in _3p_to_5p: _3p_to_5p[e2_s] = []
                if not e1_e in _5p_to_gene_info: _5p_to_gene_info[e1_e] = []
                if not e2_s in _3p_to_gene_info: _3p_to_gene_info[e2_s] = []
                    
                if not e2_s in _5p_to_3p[e1_e]:
                    _5p_to_3p[e1_e].append(e2_s)
                    _3p_to_5p[e2_s].append(e1_e)
                
                if not gene_inf in _5p_to_gene_info[e1_e]: _5p_to_gene_info[e1_e].append(gene_inf)
                if not gene_inf in _3p_to_gene_info[e2_s]: _3p_to_gene_info[e2_s].append(gene_inf)
                
                ss_G.add_node(e1_s, type=3)
                ss_G.add_node(e1_e, type=5)
                ss_G.add_node(e2_s, type=3)
                ss_G.add_node(e2_e, type=5)
                ss_G.add_edge(e1_e , e2_s, type='intron')
                ss_G.add_edge(e2_s , e2_e, type='exon')
                ss_G.add_edge(e1_s , e1_e, type='exon')

            else:
                #REV
                if not e2_s in _5p_to_3p: _5p_to_3p[e2_s] = []
                if not e1_e in _3p_to_5p: _3p_to_5p[e1_e] = []
                if not e2_s in _5p_to_gene_info: _5p_to_gene_info[e2_s] = []
                if not e1_e in _3p_to_gene_info: _3p_to_gene_info[e1_e] = []

                if not e1_e in _5p_to_3p[e2_s]:
                    _5p_to_3p[e2_s].append(e1_e)
                    _3p_to_5p[e1_e].append(e2_s)
                
                if not gene_inf in _5p_to_gene_info[e2_s]: _5p_to_gene_info[e2_s].append(gene_inf)
                if not gene_inf in _3p_to_gene_info[e1_e]: _3p_to_gene_info[e1_e].append(gene_inf)
                
                ss_G.add_node(e2_e, type=3)
                ss_G.add_node(e2_s, type=5)
                ss_G.add_node(e1_e, type=3)
                ss_G.add_node(e1_s, type=5) #note, last one will be not actually a 5'...
                ss_G.add_edge(e2_s , e1_e, type='intron')
                ss_G.add_edge(e1_e , e1_s, type='exon')
                ss_G.add_edge(e2_e , e2_s, type='exon')
                
        """
        add the last exon
        """
        e1_s, e1_e = self.exons[-1]
        #add the exons
        if not e1_s in _exon_s_to_e: _exon_s_to_e[e1_s] = []
        if not e1_e in _exon_e_to_s: _exon_e_to_s[e1_e] = []
        
        e1_tup = tuple([e1_s, e1_e])
        if not e1_tup in all_uniq_exons: all_uniq_exons[e1_tup] = 1
        
        if not e1_e in _exon_s_to_e[e1_s]:
            _exon_s_to_e[e1_s].append(e1_e)
        if not e1_s in _exon_e_to_s[e1_e]:
            _exon_e_to_s[e1_e].append(e1_s)
        
        """
        add gene info to the END of the last exon, 
        even though not techinically a 5'
        """
        if self.strand == 1:
            if not e1_e in _5p_to_gene_info: _5p_to_gene_info[e1_e] = []
            if not gene_inf in _5p_to_gene_info[e1_e]: _5p_to_gene_info[e1_e].append(gene_inf)
        else: 
            if not e1_s in _5p_to_gene_info: _5p_to_gene_info[e1_s] = []
            if not gene_inf in _5p_to_gene_info[e1_s]: _5p_to_gene_info[e1_s].append(gene_inf)
        
    @classmethod
    def init_from_feature(cls, contig, feature):
        """
        assuming that all exons are orders s<e
        """
        gene_name = "gene_name" in feature.qualifiers and feature.qualifiers['gene_name'][0] or "."
        kwargs = {'contig' : contig,
                  'feature_ID': feature.id,
                  'gene_name': gene_name,
                  'gene_ID': feature.qualifiers['geneID'][0],
                  'strand': feature.strand,
                  'g_start': feature.location.start.position,
                  'g_end': feature.location.end.position,
                  'exons':[[s.location.start.position,
                            s.location.end.position] for s in feature.sub_features if s.type=="exon"]
                  }
        return cls(**kwargs)


def fetch_trascripts(GFF_parser, limit_info={}):
    """
    NOT IMPLEMENTED
    limit_info={"gff_id": ["13"]}
    """
    exit(1)
    GFF_parser = GFF.GFFParser()
    for rec in GFF_parser.parse_in_parts(open(o.fn_gff), limit_info=limit_info):
    #iterate over chromosomes
        contig_seq = fa.get_sequence(rec.id)
        
        D_to_A, A_to_D = get_donor_acceptor_pairs(rec)
        for feature in rec.features:
            if feature.type in ["mRNA", "transcript"]:
                t = transcript.init_from_feature(feature)
                alt_ss = t.get_all_3pSS(contig_seq)
                print alt_ss
