from __future__ import print_function

import itertools
import numpy as np
import pandas as pd
import sys
from collections import defaultdict  
from scipy.stats import fisher_exact 
from scipy.stats import chi2_contingency
import pdb
import time

from Queue import Queue
import multiprocessing as mp
import math

import RNA

class fold_process(mp.Process):

    def __init__(self,):
        mp.Process.__init__(self)
        self.exit = mp.Event()

    def terminate(self):
        print("{proc} terminating".format(proc=self.p_num), file=sys.stderr)
        self.exit.set()
   
    def set_params(self, p_num, q_input, q_output, n_jobs):
        
        self.p_num = p_num
        self.q_input = q_input
        self.q_output = q_output
        self.n_jobs = n_jobs
        
    def run(self):

        n_done = 0
        folds = []
        while n_done < self.n_jobs:
            seq_dict = self.q_input.get()
            if seq_dict == None:
                assert False, "n_folded:%d n_folds:%d"%(n_folded, self.n_folds)
                break
            
            fold, e = RNA.fold(seq_dict['seq'])
            seq = seq_dict['seq']
            gc = (seq.count("G") + seq.count("C"))/float(len(seq))

            d = {"fold":fold,
                 "e":e,
                 "gc":gc}

            d.update(seq_dict)
            folds.append(d)
            n_done +=1
            self.q_input.task_done()
        
        self.q_output.put(folds) 
        self.terminate()

class folder():

    def __init__(self):
        pass

    def do_fold(self, seq_dicts, n_procs = 64):
        """
        pass in list of seqs
        return folds
        """
        
        assert len(seq_dicts) > n_procs

        q_input = mp.JoinableQueue()
        q_output = mp.JoinableQueue()
        
        n_procs_remain = n_procs
        n_jobs_total = len(seq_dicts)
        n_jobs_remain = n_jobs_total
        procs = [] 
        
        for i in range(n_procs):
            n_proc_jobs = int(math.ceil(n_jobs_remain/float(n_procs_remain)))
            n_procs_remain-=1
            n_jobs_remain-=n_proc_jobs
            
            p=fold_process()
            p.set_params(i, q_input, q_output, n_proc_jobs)
            p.daemon=True
            p.start()
            print("process %d started"%i, file=sys.stderr)
            procs.append(p)

        curr_count = 0
        for seq_dict in seq_dicts:
            q_input.put(seq_dict)
        
        while not q_input.empty():
            time.sleep(10)
            print(q_input.qsize(), n_jobs_total, file=sys.stderr)
        
        self.folds = []

        for p in procs:
            print("extracting from queue... ", file=sys.stderr)
            _folds = q_output.get()
            q_output.task_done()
            self.folds.extend(_folds)
            print("done", file=sys.stderr)

