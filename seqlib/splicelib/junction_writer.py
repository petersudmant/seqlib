



class JunctionWriter(object):
    """
    JUNCTIONS are the last base and first base of the left and right exons
    respectively, HOWEVER, STAR and tophat request 0 based and 1 based 
    coordinates respectively:

    tophat:
        http://ccb.jhu.edu/software/tophat/manual.shtml
        <chrom> <left> <right> <+/->
        left and right are zero-based coordinates, and specify the last 
        character of the left sequenced to be spliced to the first character 
        of the right sequence, inclusive. That is, the last and the first 
        positions of the flanking exons. Users can convert junctions.bed 
        (one of the TopHat outputs) to this format using 
        bed_to_juncs < junctions.bed > new_list.juncs 
        where bed_to_juncs can be found under the same folder as tophat
    
    rnaSTAR:
        https://code.google.com/p/rna-star/downloads/detail?name=STARmanual_2.3.0.1.pdf
        STAR can also utilize annotations formatted as a list of splice junctions 
        coordinates in a text file: --sjdbFileChrStartEnd /path/to/sjdbFile.txt. 
        This file should contains 4 columns separated by tabs:
            Chr \tab Start \tab End \tab Strand=+/-/.
        Here Start and End are first and last bases of the introns 
        (1-based chromosome coordinates). This file can be used in addition to the 
        --sjdbGTFfile, in which case STAR will extract junctions from both files.
        Note, that the --sjdbFileChrStartEnd file can contain duplicate (identical) 
        junctions, STAR will collapse (remove) duplicate junctions.

    THUS we output X.tophat.juncs and X.STAR.juncs
    
    0 based - 012345678
    1 based - 123456789       
    genome  - XGT---AGX

    tophat chrXX\t0\t8\t+
    STAR chrXX\t2\t8\t+
    
    """
    def __init__(self, fn_prefix):
        
        fn_tophat="%s.tophat.juncs"%fn_prefix
        fn_STAR="%s.tophat.juncs"%fn_prefix
        self.F_tophat = open(fn_prefix,'w')
        self.F_STAR = open(fn_STAR,'w')
        
        self.junc_tups_output = {}

    def write(self, junc_tups):

        o_tophat_junc_lines = []
        o_STAR_junc_lines = []
        pattern = "{contig}\t{left}\t{right}\t{strand}"
        
        for tup in junc_tups:
            if not tup in self.junc_tups_output:
                o_tophat_junc_lines.append(pattern.format(contig=tup[0],
                                                          left = tup[1],
                                                          right= tup[2],
                                                          strand = tup[3]))
                o_STAR_junc_lines.append(pattern.format(contig=tup[0],
                                                        left = tup[1]+2,
                                                        right= tup[2],
                                                        strand = tup[3]))
                self.junc_tups_output[tup] = True
        self.F_tophat.write("%s\n"%("\n".join(o_tophat_junc_lines)))
        self.F_STAR.write("%s\n"%("\n".join(o_STAR_junc_lines)))

