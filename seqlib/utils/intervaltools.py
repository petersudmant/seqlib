from bx.intervals.intersection import Interval, IntervalTree
from bx.intervals.cluster import ClusterTree

class ClusterTool(object):
    """
    get overlapping / non-overlapping clusters
    """
    def __init__(self):
        self.idx = 0
        self.objects_by_idx = {}
        self.cluster_trees_by_contig = {}

    def get_next_idx(self):
        self.idx+=1
        return self.idx
    
    def insert(self, contig, s, e, item=None):
        idx = self.get_next_idx()
        self.objects_by_idx[idx] = item
        if not contig in self.cluster_trees_by_contig:
            self.cluster_trees_by_contig[contig] = ClusterTree(0,0) 
        
        self.cluster_trees_by_contig[contig].insert(s,e,idx)
                    
    def get_regions(self):
        ret = {}
        for contig, tree in self.cluster_trees_by_contig.items():
            ret_list = []
            for cluster_tup in tree.getregions():
                s,e, idxs = cluster_tup
                ret_list.append((s,e,[self.objects_by_idx[idx] for idx in idxs]))

            ret[contig] = ret_list
        
        return ret

class IntervalTool(object):
    """
    check if something overlaps a region quickly
    """

    def __init__(self):
        self.intervals_by_contig = {}

    def insert(self, contig, s, e, item=None):

        if not contig in self.intervals_by_contig:
            self.intervals_by_contig[contig] = IntervalTree()
        
        self.intervals_by_contig[contig].insert_interval(Interval(s,e, item))
    
