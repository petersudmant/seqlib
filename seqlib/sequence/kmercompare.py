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


class count_kmers_process(mp.Process):

    def __init__(self,):
        mp.Process.__init__(self)
        self.exit = mp.Event()

    def terminate(self):
        print("{proc} terminating".format(proc=self.p_num), file=sys.stderr)
        self.exit.set()
   
    def set_params(self, p_num, q_input, q_output, seq1, seq2, l1, l2, n_jobs):
        self.p_num = p_num
        self.q_input = q_input
        self.q_output = q_output 
        self.s1 = seq1
        self.s2 = seq2
        self.l1 = l1
        self.l2 = l2
        self.n_jobs = n_jobs
        
    def run(self):

        n_done = 0
        kmer_counts = []
        while n_done < self.n_jobs:
            kmer = self.q_input.get()
            if kmer == None:
                assert False, "n_folded:%d n_folds:%d"%(n_folded, self.n_folds)
                break
            
            c1 = self.s1.count(kmer)
            c2 = self.s2.count(kmer)
            #OR, p_f = fisher_exact([[c1,self.l1], [c2,self.l2]])
            chi, p_chi, dof, exp = chi2_contingency([[c1,self.l1], [c2,self.l2]])
            if c2==0: c2=1 
            OR = (c1/float(self.l1))/(c2/float(self.l2))
            kmer_counts.append({"s1":c1,
                                "s2":c2,
                                "t1":self.l1,
                                "t2":self.l2,
                                "OR":OR,
                                "p":p_chi, 
                                "kmer":kmer})
            n_done +=1
            self.q_input.task_done()
        
        self.q_output.put(kmer_counts) 
        self.terminate()

class KmerCompare():

    def __init__(self, kmer_len=6):
        self.kmer_gen = itertools.product("ATCG",repeat=kmer_len)
        self.kmer_len = kmer_len

    def count_kmers(self, seqs1, seqs2, n_procs = 16):
        """
        pass in two lists of stringsa
        return counts
        """
        
        self.seqs1 = seqs1
        self.seqs2 = seqs2
        
        l1 = np.sum([len(s)-self.kmer_len for s in seqs1])
        l2 = np.sum([len(s)-self.kmer_len for s in seqs2])
        
        s1 = "|".join(seqs1).upper()
        s2 = "|".join(seqs2).upper()
    
        q_input = mp.JoinableQueue()
        q_output = mp.JoinableQueue()
        
        n_procs_remain = n_procs
        n_jobs_total = 4**self.kmer_len
        n_jobs_remain = n_jobs_total
        procs = [] 
        
        for i in range(n_procs):
            n_proc_jobs = int(math.ceil(n_jobs_remain/float(n_procs_remain)))
            n_procs_remain-=1
            n_jobs_remain-=n_proc_jobs
            
            p=count_kmers_process()
            p.set_params(i, q_input, q_output, s1, s2, l1, l2, n_proc_jobs)
            p.daemon=True
            p.start()
            print("process %d started"%i, file=sys.stderr)
            procs.append(p)

        curr_count = 0
        for kmer in self.kmer_gen:
            kmer_str = "".join(kmer)
            q_input.put(kmer_str)
        
        while not q_input.empty():
            time.sleep(10)
            print(q_input.qsize(), n_jobs_total, file=sys.stderr)
        
        self.kmer_counts = []

        for p in procs:
            print("extracting from queue... ", file=sys.stderr)
            _kmer_counts = q_output.get()
            q_output.task_done()
            self.kmer_counts.extend(_kmer_counts)
            print("done", file=sys.stderr)

    def get_kmer_stats(self, kmer_count):
        kmer = kmer_count['kmer'] 
        s1_cs = np.array([seq.count(kmer) for seq in self.seqs1])
        s2_cs = np.array([seq.count(kmer) for seq in self.seqs2])
        
        mu_s1, med_s1, mx_s1 = np.mean(s1_cs), np.median(s1_cs), np.amax(s1_cs)
        mu_s2, med_s2, mx_s2 = np.mean(s2_cs), np.median(s2_cs), np.amax(s1_cs)
        n_s1_gt0 = np.sum(s1_cs>0) 
        n_s2_gt0 = np.sum(s2_cs>0) 
        
        return {"mu_s1_count":mu_s1,
                "med_s1_count":med_s1,
                "max_s1_count":mx_s1,
                "n_s1_w_kmer":n_s1_gt0,
                "mu_s2_count":mu_s2,
                "med_s2_count":med_s2,
                "max_s2_count":mx_s2,
                "n_s2_w_kmer":n_s2_gt0}



    def get_significant_kmers(self, fdr = 0.01, Yekuteili=True):
        self.kmer_counts  = sorted(self.kmer_counts, key = lambda x: x['p'])
        m = float(len(self.kmer_counts))
        cm=1
        if Yekuteili:
            """
            correct for non-independence of tests 
            """
            cm =  np.sum(1/(np.arange(m)+1))

        outrows = [] 
        for i, kmer_count in enumerate(self.kmer_counts):
            p = kmer_count['p']
            j=i+1
            if p<=((j/(m*cm))*fdr):
                BHY_sig = True
            else:
                BHY_sig = False

            kmer_stats = self.get_kmer_stats(kmer_count)
            kmer_count.update(kmer_stats)
            kmer_count.update({"BHY-significant":BHY_sig})
            outrows.append(kmer_count)
        
        return pd.DataFrame(outrows)



        
