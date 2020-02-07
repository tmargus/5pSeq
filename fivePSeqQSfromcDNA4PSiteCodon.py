#!/usr/bin/env python
# Python 3
#
# Script for calculating Queuing Score in a middle of mRNA starting from specified codon in P-Site
#
# WARNIGS! dont use it as it is !!!!!
#

__author__		= "Tonu Margus"
__copyright__	= "Copyright 2020"
__version__		= "0.0.1"
__email__		= "tonu.margus@gmail.com"
__status__		= "beta"

import re
import sys
import pysam
import argparse
import pandas as pd
import pickle
# 
# -i input comes from the script fivePSeqMap2cdnawitfOffset.py (https://github.com/tmargus/5pSeq)
#
parser = argparse.ArgumentParser(description="P-site assigned reads (offset 15)" )
parser.add_argument('-i',  type=str, help='input py3 pickle file from "fivePSeqMap2cdnawithOffset.py" ')
parser.add_argument('-f',  type=str, help='')
parser.add_argument('-o',  type=str, help='output_QS_.csv offset added lik QS_./QS_offset.', default='output_QS_.csv')
args = parser.parse_args()

## set some filters
min_counts   = 10  # minimal number of reads mapped to gene
min_aver_cov = 0.2 # minimal average coverage  all_counts/gene_length_in_codons
 ##
sys.stderr.write("\n\
-i  input sorted BAM : {}\n\
-f  cDNA in fasta    : {}\n\
-o  output pkl       : {}\n\n".format(args.i, args.f, args.o))

usage = 'python fivePSeqQSfromcDNA4PSiteCodon.py  -i experiment_cond.pkl -o  output_cdna.csv"'

if args.i==None:
     sys.exit("\n  usage:\n\t{}\n".format(usage))
###################################################################



###################################################################
# Functions
# (4) QS for coding region
def codon_positions(df, codon):
    '''return list of positions for codon
    '''
    n  = df.shape[0]
    df = df.loc[60:n-30,]
    return list(df.loc[df['codon']==codon,].index)

def qsW(df, i, per=10):
    ''' weighted sum of signal
    df: dataframe
    i: position for codon/amino acid in P-Site
    per: period for wave in codons/amino acids
    '''
    s1 = df.loc[i+1,'counts']          # mono-
    s2 = df.loc[i-per+1,'counts']  *2  # di-
    s3 = df.loc[i-per*2+1,'counts']*3  # tri-
    
    return s1+s2+s3

def qs_bkg(df,i,per=10):
    ''' background between given period
    df: dataframe
    i: position for codon/amino acid in P-Site
    per: period for wave in codons/amino acids
    '''
    
    b1 = df.loc[i-per+1:i-1,'counts'].mean()         # mono- and disomes
    b2 = df.loc[i-per*2+1:i-per-1,'counts'].mean()   # di- and trisomes
    b3 = df.loc[i-per*3+1:i-per*2-1,'counts'].mean() # tri- and tetrasomes
    
    return (b1+b2+b3)/3

def block(index_l, i, incr=1):
    ''' returns number > 1 for block 
    of identical elements lined up
    '''
    a = np.array(index_l)
    return len(a[(a>=i-incr) &(a<=i+incr)])


###################################################################
# main program
# 
f = "WT20C_V_QS-coding-region_AGG_P-site.csv" ############outfile 
of = open(f, "w")
l = []

report = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format("pos","QS","qs_sum","bkg","E","P","A","rep","cDNA")
of.write(report +"\n")

c = c1 = 0

for k in list(d.keys()): 

    df   = d[k]
    cp_l = codon_positions(df, 'AGG')
    
    for i in cp_l:
        qs_sum = qsW(df, i, per=10)
        bkg    = qs_bkg(df,i, per=10)
        if bkg < 0.05:   # avoid div 0
            c1+=1
            continue
            
        QS     = qs_sum/bkg
        A      = df.loc[i+1, 'codon']
        P      = df.loc[i, 'codon']
        E      = df.loc[i-1, 'codon']
        comment= '*' * block(cp_l, i, incr=1)
        report = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(i,QS,qs_sum,bkg,E,P,A,comment,k)
        of.write(report +"\n")
        #print("{:7,}  {:6.2f}  {:6.2f}  {:6.2f}  {:4} {:4} {:4} {:4} {:16}".format(i,QS,qs_sum,bkg,E,P,A,comment,k))
        c+=1
        
of.close()
print(f)
