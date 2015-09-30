from pysam.calignmentfile cimport AlignmentFile, AlignedSegment
from pysam.calignmentfile cimport IteratorColumn
from pysam.calignmentfile cimport IteratorColumnRegion
from pysam.calignmentfile cimport IteratorRow
from pysam.calignmentfile cimport PileupColumn, PileupRead
from pysam.ctabix cimport Tabixfile
cimport numpy as np
import numpy as np

cdef AlignmentFile samfile
cdef Tabixfile tabixfile

#@cython.boundscheck(False)
#@cython.wraparound(False)

def get_seg_cvg(IteratorColumnRegion col_iter, int shape):
    
    cdef IteratorRow read_iter
    cdef int n_total = 0 
    cdef PileupColumn col
    cdef PileupRead read
    
    cdef Py_ssize_t i = 0
    cdef np.ndarray cvg = np.zeros(shape, dtype=np.int)
    
    for col in col_iter:
        n_total = 0
        for read in col.pileups:
            if  read.is_del==0:
                n_total+=1
        cvg[i]=n_total
        i+=1
    
    return cvg

def testCountBAM(AlignmentFile samfile):
    '''test reading from a BAM file accessing
    the flag field directly.'''

    cdef AlignedSegment read
    cdef int n = 0
    
    for read in samfile.fetch():
        flag = read._delegate.core.flag
        n += 1
            
    return n

def testCountGTF(Tabixfile tabixfile):
    '''test reading from a tabixfile.'''
    
    cdef int n = 0

    for entry in tabixfile.fetch():
        n += 1

    return n
