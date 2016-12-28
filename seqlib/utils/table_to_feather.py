import argparse
import pandas as pd
import feather


if __name__=="__main__":
    parser = argparse.ArgumentParser()    
    parser.add_argument("--fn_input")
    parser.add_argument("--fn_output")
    parser.add_argument("--compression", default=None)
    parser.add_argument("--column_names", default=None, n_args="+")
    o = parser.parse_args()
    
    if o.column_names==None:
        t = pd.read_csv(o.fn_input, header=0, sep="\t", compression=o.compression)
    else:
        t = pd.read_csv(o.fn_input, 
                        header=0, 
                        sep="\t", 
                        compression=o.compression, 
                        header=None,
                        names=o.column_names)
    
    feather.write_dataframe(t, o.fn_output)
