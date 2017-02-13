import argparse
import pandas as pd
import feather
import pdb

if __name__=="__main__":
    parser = argparse.ArgumentParser()    
    parser.add_argument("--fn_input", nargs="+")
    parser.add_argument("--fn_output")
    parser.add_argument("--compression", default=None)
    parser.add_argument("--column_names", default=None, nargs="+")
    o = parser.parse_args()
    
    tables = []
    if o.column_names==None:
        for fn in o.fn_input:
            t = pd.read_csv(fn, header=0, sep="\t", compression=o.compression)
            tables.append(t)
    else:
        for fn in o.fn_input:
            t = pd.read_csv(fn,
                            sep="\t", 
                            compression=o.compression, 
                            header=None,
                            names=o.column_names)
            tables.append(t)
    t = pd.concat(tables) 
    feather.write_dataframe(t, o.fn_output)
