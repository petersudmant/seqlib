
class JunctionCount(object):
    """
    if readlen == 5 
    ____|   |____
    4567|   |
      ---    --
    
    """
    def __init__(self, junc, readlen):
        self.junc_l = junc[1]
        self.junc_r = junc[2]
        
        self.wnd_start = self.junc_l - (readlen-2)
        self.readlen = readlen
        self.read_counts = np.zeros(readlen-1)
    
    def add_read(self, read_s):
        self.read_counts[read_s-self.wnd_start]+=1
    
    def calc_entropy(self, min_o=1):
        """
        min_overlap MUST be at least 1
        """
        assert min_o>=1
        l = self.read_counts.shape[0]
        c_read_counts = self.read_counts[min_o-1:l+1-min_o]
        t = np.sum(c_read_counts)
        if t==0: return 0
        ps = c_read_counts / t
        lps = np.log(ps)
        lps[ps==0]=0
        e = -np.sum((ps*lps)/np.log(2))
        return e
    
    def n_supporting_reads(self, min_o=1): 
        """
        min_overlap MUST be at least 1
        """
        assert min_o>=1
        l = self.read_counts.shape[0]
        c_read_counts = self.read_counts[min_o-1:l+1-min_o]
        t = np.sum(c_read_counts)
        return t


    def min_lr_overhang(self):
        """
        if readlen == 5 
        ____|   |____
        4567|   |
          ---    --
        0123
           -----
        """
        w_reads = np.where(self.read_counts!=0)[0]
        if w_reads.shape[0]==0: return 0
        max_l = self.readlen-np.amin(w_reads)+1
        max_r = np.amax(w_reads)+1
        return min(max_l, max_r)
