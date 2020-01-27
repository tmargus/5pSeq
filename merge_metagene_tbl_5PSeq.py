#!/usr/bin/env python
#
# merges metagene tables computed with   `hdf_2_metagene_tables_stop_5PSeq.py` ->  single table
#

__author__		= "Tonu Margus"
__copyright__	= "Copyright 2019"
__version__		= "0.0.1"
__email__		= "tonu.margus@gmail.com"
__status__		= "beta"

import re
import sys
import glob
import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='Merges metagene tables of different samples')
parser.add_argument('-files', type=str, help='single sample waviness tables to merge: tray "ls -1 *.csv"')
parser.add_argument('-o', type=str, help='output file name - defualt is "output_tbl.csv"', default="output_tbl.csv")
parser.add_argument('-col',  type=str, help='column name for coverage: "raw"; "rpm"; "normalised"(def)', default="normalised")
args = parser.parse_args()


print("\n\
-files   : {}\n\
-o       : {}\n".format(args.files, args.o))

usage = ' ./merge_metagene_tbl.py  -files "8-MetagTbl/*sum.csv" -o metagene_summ_tbl.csv'


if (args.files==None)|(args.o==None):
     sys.exit("\n  usage:\n\t{}\n".format(usage))

flist   = glob.glob(args.files)
flist.sort()
outfile = args.o
col     = args.col

def read_series_in(file, column, name):
    """ returns pandas series
    """
    df = pd.read_csv(f, sep='\t', index_col=0)
    s  = df[column]
    s.name = name
    return s

dfSum = pd.DataFrame()

for f in flist:
    #sys.stderr.write("{}\n".format(f))
    m  = re.split('/|-',f)           #### might vary from input
    n  = "{}-{}".format(m[2],m[3])   #### might vary from input
    s  =  read_series_in(f, col, n)
    
    if dfSum.empty:
        dfSum = pd.DataFrame(s, index=s.index)
    else:
        dfSum[s.name] = s

dfSum.to_csv(outfile, sep="\t")
print("df.shape:  {}".format(dfSum.shape))
print("Columns:   {}".format(list(dfSum.columns)))
print(outfile)
dfSum[:3] 