



class JunctionWriter(object):
    """
    output junctions in two forms
        STANDARD X.juncs
            STANDARD is bed type coordinates, 0-based, half open
        1_based X.1_based.juncs
            1_based is the 1 based coords of the first bp of the intron and 
            the last bp of the intron for STAR
    """
    def __init__(self, fn):
        
        fn_1_based = fn.replace("juncs","1_based.juncs")
        self.junc_tups_output = {}
        self.F = open(fn, 'w')
        self.F_1_based = open(fn_1_based, 'w')
    
    def write(self, junc_tups):

        o_junc_lines = []
        o_junc_1_based_lines = []
        pattern = "{contig}\t{left}\t{right}\t{strand}"
        
        for tup in junc_tups:
            if not tup in self.junc_tups_output:
                o_junc_lines.append(pattern.format(contig=tup[0],
                                                   left = tup[1],
                                                   right= tup[2],
                                                   strand = tup[3]))
                #one based have the same right coordinate, but a left coordinate +1
                o_junc_1_based_lines.append(pattern.format(contig=tup[0],
                                                   left = tup[1]+1,
                                                   right= tup[2],
                                                   strand = tup[3]))
                self.junc_tups_output[tup] = True
        self.F.write("%s\n"%("\n".join(o_junc_lines)))
        self.F_1_based.write("%s\n"%("\n".join(o_junc_1_based_lines)))

