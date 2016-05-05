import argparse
import pandas as pd
import feather


if __name__=="__main__":
    parser = argparse.ArgumentParser()    
    parser.add_argument("--fn_input")
    parser.add_argument("--fn_output")
    parser.add_argument("--compression", default=None)
    o = parser.parse_args()
    
    t = pd.read_csv(o.fn_input, header=0, sep="\t", compression=o.compression)
    feather.write_dataframe(t, o.fn_output)
