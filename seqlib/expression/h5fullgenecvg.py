import tables
from coveragedata import CoverageData 
import numpy as np


class h5FullGeneCvg(object):
    def __init__(self, fn):
        self.h5 = tables.openFile(fn, mode='w')
        
        filt = tables.Filters(complevel=5, complib='blosc')
        self.a_cvg = self.h5.createEArray(self.h5.root, 
                                          'cvg', 
                                          tables.Float32Atom(), 
                                          shape=(0,), 
                                          filters=filt,
                                          expectedrows=65000000)
        
        self.a_idx = self.h5.createEArray(self.h5.root, 
                                          'idx', 
                                          tables.UInt16Atom(), 
                                          shape=(0,), 
                                          filters=filt,
                                          expectedrows=65000000)
        
        self.a_GeneID = self.h5.createEArray(self.h5.root, 
                                             'geneID', 
                                             tables.StringAtom(20),
                                             shape=(0,), 
                                             filters=filt,
                                             expectedrows=20000)
        
        self.a_start = self.h5.createEArray(self.h5.root, 
                                            'start', 
                                            tables.UInt16Atom(),
                                            shape=(0,), 
                                            filters=filt,
                                            expectedrows=20000)
        
        self.a_stop = self.h5.createEArray(self.h5.root, 
                                           'stop', 
                                           tables.UInt16Atom(),
                                           shape=(0,), 
                                           filters=filt,
                                           expectedrows=20000)
        
        self.a_length = self.h5.createEArray(self.h5.root, 
                                             'length', 
                                             tables.UInt16Atom(),
                                             shape=(0,), 
                                             filters=filt,
                                             expectedrows=20000)
        self.curr_idx = 0
        self.GeneIDs = []
        self.lengths = []
        self.starts = []
        self.stops = []
    
    def extend(self, cvg_obj):

        self.a_cvg.append(cvg_obj.RNA_cvg_view)
        self.a_idx.append(np.ones(cvg_obj.RNA_cvg_view.shape[0])*self.curr_idx)
        
        self.GeneIDs.append(cvg_obj.g.gene_id)
        self.lengths.append(cvg_obj.RNA_cvg_view.shape[0])
        self.starts.append(cvg_obj.RNA_cvg_view_START)
        self.stops.append(cvg_obj.RNA_cvg_view_STOP)

        self.curr_idx+=1
    
    def close(self):
        
        self.a_GeneID.append(np.array(self.GeneIDs))
        self.a_length.append(np.array(self.lengths))
        self.a_start.append(np.array(self.starts))
        self.a_stop.append(np.array(self.stops))

        self.h5.close()

