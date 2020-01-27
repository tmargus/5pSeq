#!/usr/bin/env python
# 
# Pools two replicates together 5PSeq
# input tables in HDF  format  
# doubleindexed "For_rpm/I"  "Rev_rpm/I"
# output in HDF format


__author__		= "Tonu Margus"
__copyright__	= "Copyright 2020"
__version__		= "0.1"
__email__		= "tonu.margus@gmail.com"
__status__		= "beta"

import re
import sys
import glob
import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='Merges metagene tables of different samples')
parser.add_argument('-i1', type=str, help='first sample S1_idx_iv.h5')
parser.add_argument('-i2', type=str, help='second sample S2_idx_iv.h5')
parser.add_argument('-o', type=str, help='output pooled S1S2_idx_iv.h5')
parser.add_argument('-col',  type=str, help='column name for coverage: "raw"; "rpm"; "normalised"(def)', default="rpm")
args = parser.parse_args()

print("\n\
-i1   : {}\n\
-i2   : {}\n\
-col  : {}\n\
-o    : {}\n".format(args.i1, args.i2, args.col, args.o))

usage = ' ./pool_replicates_h5.py -i1  6-AssignRaw/S1_idx_iv.h5 -i2 6-AssignRaw/S1_idx_iv.h5 -o S1S2_idx_iv_pooled.h5'


# test 
if (args.i1==None)|(args.i2==None)|(args.o==None):
     sys.exit("\n  usage:\n\t{}\n".format(usage))

# open files
try:
    hd5_1 = pd.HDFStore(args.i1, "r")
except:
    print('Could not open input "{}"'.format(args.i1))
    sys.exit(1) 

try:
    hd5_2 = pd.HDFStore(args.i2, "r") 
except:
    print('Could not open input "{}"'.format(args.i2))
    sys.exit(1) 


store = pd.HDFStore(args.o, complevel=5, complib="zlib", mode="w")

col= args.col

for k in hd5_1.keys():
    #sys.stderr.write("{:10}...\n".format(k))
    columns = hd5_1[k].columns
    Chr     = hd5_1[k].iloc[0]['Chr']
    Strand  = hd5_1[k].iloc[0]['Strand']
    # get values and sum
    s1 = hd5_1[k][col]
    s2 = hd5_2[k][col]
    s  = s1.add(s2, fill_value=0)
    s.name=col
    # make df
    
    df = pd.DataFrame(s, index=s.index)
    df['Chr'] = Chr
    df['Strand']= Strand
    df = df[columns]
    
    store.put(k, df)

hd5_1.close()
hd5_2.close()
store.close()
print("//")
#print("\n OutFile: {}\n".format(args.o))
