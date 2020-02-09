#!/usr/bin/env python
# Python 3
#
# Script for calculating Queuing Score in a middle of mRNA starting from specified codon in P-Site
#
# WARNIGS! dont use it as it is !!!!!
#

__author__		= "Tonu Margus"
__copyright__	= "Copyright 2020"
__version__		= "0.0.2"
__email__		= "tonu.margus@gmail.com"
__status__		= "beta"

import re
import sys
import argparse
import pandas as pd
import numpy as np
import pickle
# 
# -i input comes from the script fivePSeqMap2cdnawitfOffset.py (https://github.com/tmargus/5pSeq)
#
parser = argparse.ArgumentParser(description="Queuing scores for fixed P-site {61} A-Site pairs in the middle of coding sequence" )
parser.add_argument('-i',  type=str, help='input py3 pickle file from "fivePSeqMap2cdnawithOffset.py" ')
parser.add_argument('-codon',  type=str, help='P-Site codon', default="AAA")
parser.add_argument('-per',  type=int, help='Periodicity in codons for QS; default is 10 codons', default=10)
parser.add_argument('-o',  type=str, help='output_QS.csv', default='QS-coding-region_AAA_P-site.csv')
args = parser.parse_args()

## set some filters
#min_counts = 10  # minimal number of reads mapped to gene
min_bkg    = 0.05 # minimal average coverage  all_counts/gene_length_in_codons
 ##
sys.stderr.write("\n\
-i     input pickled dict : {}\n\
-codon     P-Site codon   : {}\n\
-per periodicity in codons: {}\n\
-o     output_table.csv   : {}\n\n".format(args.i, args.codon, args.per, args.o))

usage = 'python fivePSeqQSfromcDNA4PSiteCodon.py  -codon AAA -i experiment_cond.pkl -o  output_cdna.csv"' 

if args.i==None:
     sys.exit("\n  usage:\n\t{}\n".format(usage))
     
###################################################################

codon = args.codon
of    = open( args.o, "w")
l     = []

report = "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format("pos","QS","qs_sum","bkg","E","P","A","rep","cDNA")
of.write(report +"\n")

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
    s4 = df.loc[i-per*3+1,'counts']*4  # tri-
    
    return s1+s2+s3+s4

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

with open(args.i, 'rb') as handle:
    d = pickle.load(handle)
       
print("{:,} cDNAs from input file".format(len(d)))
print("QS for P-Site {}  A-Site NNN pairs".format(codon)) 

c = c1 = cc = 0

for k in list(d.keys()): 

    df   = d[k]
    cp_l = codon_positions(df, codon=codon)
    
    for i in cp_l:
        if c == 10000:
            cc+=c
            print("{:,} ...".format(cc))
            c=0
            
        qs_sum = qsW(df, i, per=10)
        bkg    = qs_bkg(df,i, per=10)
        if bkg < min_bkg:   # avoid div 0
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
print("{:,} total QS for codon pairs {}[NNN]!".format(cc+c, codon))
print(args.o)
