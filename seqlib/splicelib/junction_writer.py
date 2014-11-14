



class JunctionWriter(object):
    
    def __init__(self, fn):
        
        self.junc_tups_output = {}
        self.F = open(fn, 'w')
    
    def write(self, junc_tups):

        o_junc_lines = []
        pattern = "{contig}\t{left}\t{right}\t{strand}"
        
        for tup in junc_tups:
            if not tup in self.junc_tups_output:
                o_junc_lines.append(pattern.format(contig=tup[0],
                                                   left = tup[1],
                                                   right= tup[2],
                                                   strand = tup[3]))
                self.junc_tups_output[tup] = True
        self.F.write("%s\n"%("\n".join(o_junc_lines)))

