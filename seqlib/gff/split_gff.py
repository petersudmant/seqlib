#!/bin/bash

import argparse
import sys
import pandas as pd
import numpy as np

# edited by SL, 2022/06/16
# @stacy-l on GitHub

# This is a short utility script that takes in a gff file that contains extra
# key-value pair info (ID, CDS ID, etc) in the optional 8th field. 
# This 8th field can contain a variable number of key-value pairs, so this
# script opts to just create a new column for each key-value pair, for the
# maximum # of key-value pairs observed in any given gff row.

# This script was originally written for .gff files where the 8th column's 
# first key-value pair is in the form of 'ID=<text>'. This script will will 
# do the above column-formatting and strip the 'ID=' string from the generated
# ID column. The remaining generated columns are left as-is.

# This script returns the entire table as a tab-delimited txt file with the
# suffix of '_exp.txt' replacing '.gff'.

def split_gff(fn):
  # takes in a gff and splits the 8th column by semicolon delimited
  # saves it as a tab delimited .txt file with the _exp suffix
  gff = pd.read_table(fn, header = None)
  expand_col = gff[8].str.split(';', expand = True) # expansion of 8th column
  split_gff = pd.concat([gff.drop(8, axis = 1), expand_col], axis = 1) # concatenate expansion
  split_gff.columns = pd.RangeIndex(split_gff.columns.size) # reset column names
  split_gff[8] = split_gff[8].str.lstrip('ID=') # strip ID prefix
  split_gff.to_csv(fn.rstrip('.gff')+'_exp.txt', sep='\t') # tab delimited

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('fn', help='gff file name', type = str)

  try:
      args = parser.parse_args()
      split_gff(args.fn)
  except:
      parser.print_help()
      sys.exit(0)