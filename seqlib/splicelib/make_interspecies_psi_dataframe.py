import argparse
import pandas as pd

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--fn_input")
    parser.add_argument("--psi_summary_dir")
    parser.add_argument("--fn_output")
    parser.add_argument("--species", nargs="*")
    o = parser.parse_args()

    all_ts = []
    for l in open(o.fn_input):
        species, contig, ss_5p, ss_3p_prox, ss_3p_dist = l.rstrip().split()
        if species not in o.species: continue
        ss_5p, ss_3p_prox, ss_3p_dist = int(ss_5p), int(ss_3p_prox), int(ss_3p_dist)
        f = "%s/%s.psis.summary"%(o.psi_summary_dir,species)
        #f = "%s/%s.psis.all"%(o.psi_summary_dir,species)
        t = pd.read_csv(f,header=0,index_col=None, delimiter="\t")
        t = t[t["ss_5p"] == ss_5p]
        t['species'] = species
        all_ts.append(t)
    
    t = pd.concat(all_ts)
    t.to_csv(o.fn_output, sep="\t", index=False)
    
