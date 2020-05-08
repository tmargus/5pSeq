#!/usr/bin/env python
#
# works fine
# Features to add:
#       a) split by stop codon and  
#
# Computing Queing score is changed. 
# assumes 5P-Seq uncorrected 5' positions what is mapped to A-Site with offset 17
# Peaks positions 1   -18:-16    *1  (mono-somes  STOP codon)
#                 2   -48:-46    *2    (di-somes)
#                 3   -78:-76    *3   (tri-somes)
#                 4  -108:-106   *4 (tetra-somes) . # ver. 0.1.6
# divide with avarage codon coverage   and  sum()  = QS
# 
# bkg region can be selected between peaks (-bkg 2) or before tetrasome peak area (-bkg 1 (default) ) # ver. 0.1.7
#
# -w weigthing tetra-,tri-, and disomes  YES/NO  default NO # ver. 0.1.8

__author__		= "Tonu Margus"
__copyright__	= "Copyright 2020"
__version__		= "0.1.8"
__email__		= "tonu.margus@gmail.com"
__status__		= "beta"

import sys
import argparse
import pandas as pd
import pysam


parser = argparse.ArgumentParser(description='Computes ribosomes queing score at stop')
parser.add_argument('-i', type=str, help='input table of gene coverage in *.hd5 format')
parser.add_argument('-o', type=str, help='output file name  *.csv')
parser.add_argument('-annot',  type=str, help='GTF annotation file', default='0-References/genome.gtf.gz')
parser.add_argument('-th1',  type=float, help='Summary gene coverage 150 nt before stop - 10(rpm) default', default=15)
parser.add_argument('-th2',  type=float, help='Background coverage - codon mean from -115 up to Span', default=0.15)
parser.add_argument('-span',  type=int, help='Positions before - stop recommended 150 or bigger', default=150)
parser.add_argument('-bkg',  type=int, help='region for bakground con be 1 or 2: 1 - before and 2 - between peaks', default=2)
parser.add_argument('-col',  type=str, help='column for values: "sum"; "rpm"; "counts"', default='rpm')
parser.add_argument('-w', type=str, help='Weigthing tetra-,tri-, and di-somes: YES/NO/BOTH', default="NO")
args = parser.parse_args()

message = "between peaks" if args.bkg==2 else "{} to {}".format(-args.span, -115)

print("\n\
-i           input *.h5: {}\n\
-o         output *.csv: {}\n\
-annot   annotation GTF: {}\n\
-col             column: {}\n\
-th1   region threshold: {}\n\
-th2      bkg threshold: {}\n\
-bkg      bkg region is: {}\n\
-w      weigthing peaks: {}\n\
-span    nt before stop: {}\n".format(args.i, args.o, args.annot, args.col, args.th1, args.th2, message, args.w, args.span))

usage = "./compute_queuing_5PSeq.py -i *.hdf  -o *.csv"

if (args.i==None)|(args.o==None):
     sys.exit("\n  usage:\n\t{}\n".format(usage))

infile_h5  = args.i
outfile    = args.o
thres_cover= args.th1 
thres_bkg  = args.th2
Span       = args.span     
annotation = args.annot
col        = args.col
bkg_region = args.bkg
weighted   = str.upper(args.w)

## check is there enough space for background
if (Span<109)  & (bkg==2):
    report = "Increase -ppan at least up to 108. Current value is {} ".format(Span)
    sys.exit("\n {} \n".format(report))
# Span - 115 must be bigger than 9 nt
if ((Span-115)<9) & (bkg==1):
     report = "{}nt is too little for backgrouns. Increase Span at least up to 124".format(Span-115)
     sys.exit("\n {} \n".format(report))

replacestr = "_th{}-{}_v4.csv".format(thres_cover, thres_bkg)
outfile = outfile.replace('.csv', replacestr)

############################
def df_framing(df1, index, columns, strand="+"):
    # create df2
    df2 = pd.DataFrame(0, index=index, columns=columns)
    df1 = df1.add(df2, fill_value=0, axis=1)
    if strand == "+":
        df1.reset_index(inplace=True)  # reset index
        return df1[columns]
    elif strand == "-":
        df1 = df1[::-1]  # reverts table
        df1.reset_index(inplace=True)  # reset index
        return df1[columns]
    else:
        # error
        print("ERROR! Expext '+'/'-' but found {} for strand".format(strand))
        
def get_part_from_gtf(annotation, reference=None, feature="CDS"):   
    tabixfile = pysam.TabixFile(annotation, parser=pysam.asGTF())
    return [gtf for gtf in tabixfile.fetch(reference=reference) if (gtf.feature == feature)]

def yeastChr():
    # Ordered yeast Chr list short names from ensembl
    return ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV',
            'XVI','Mito']

chr_length_d = { 'I':230218, 'II':813184, 'III':316620, 'IV':1531933, 'IX':439888, 'Mito':85779, 
'V':576874, 'VI':270161, 'VII':1090940, 'VIII':562643, 'X':745751, 'XI':666816, 'XII':1078177,
'XIII':924431, 'XIV':784333, 'XV':1091291, 'XVI':948066 }

def df_2_ContiniousSeries(df, index, column="rpm"):
    """ returns pd.Series with 0 values for all positions in a range of index.
    Input df is condensed, i. e. positions with values < 0 not included
    :param df:     condensed df, i. e.  don't contain rows with 0 in index
    :param index:   range of genome positions
    :param column:  column name to extract from df - default is 'sum'
    :return: Series
    """
    s1 = pd.Series(0, index=index)
    s2 = df[column]
    s3 = s1+s2
    return s3.fillna(0)

def reindex(s, name, index):
    """ Returns reindexed Series
    """
    s.name = name
    df = s.reset_index()
    df['rel_Pos'] = index
    return df.set_index('rel_Pos')[name] #returns pd.Series
#############################

c=c1=c2=c3=0# counting 
d = {}  # dictionary for collecting data
 
# input h5
hd5 = pd.HDFStore(infile_h5, "r")  
# metagene summary df
columns = hd5[hd5.keys()[0]].columns

for ref in yeastChr():
    #sys.stderr.write("{} \n".format(DEBUG! 2))
    sys.stderr.write("{:4} ...\n".format(ref))
    df_f = hd5['For_rpm/'+ ref]
    df_r = hd5['Rev_rpm/'+ ref]
    ### df 2 continious Serise
    index = list(range(0, chr_length_d[ref] + 1))
    s_f = df_2_ContiniousSeries(df_f, index, column=col)
    s_r = df_2_ContiniousSeries(df_r, index, column=col)

    stop_gtf_l = get_part_from_gtf(annotation, reference=ref, feature="stop_codon") # gtf part for stop

    for gtf in stop_gtf_l:
        #stop_codon = genome[ref][gtf.start:gtf.start+3]                      # get stop codon
        #stop_codon = stop_codon if gtf.strand=='+' else revcompl(stop_codon) # revcomp rev_strand

        #test coverage
        coverage_for = s_f.loc[gtf.start - Span:gtf.start + 3].sum()
        coverage_rev = s_r.loc[gtf.start:gtf.start + Span - 1].sum()
        coverage_rpm = coverage_for if gtf.strand == '+' else coverage_rev
        
        if coverage_rpm < thres_cover: #FILTER 1
            c1+=1
            continue
        s_1 = s_f.loc[gtf.start - Span:gtf.start + Span]      # Forvard
        s_2 = s_r.loc[gtf.end - Span - 1:gtf.end + Span - 1]  # Reverse
        # proper index & dataframe
        s = s_1 if gtf.strand == '+' else s_2[::-1]
        index = list(range(-Span, Span + 1))
        s = reindex(s, name=str(col), index=index)
        # get regions
        # codon wise 
        if bkg_region == 1:
            bkg = s.loc[-Span:-115].mean()*3      # 5PSeq & Series
        elif bkg_region == 2:
            bkg = (s.loc[-43:-21].mean()*3 + s.loc[-73:-51].mean()*3 + s.loc[-103:-81].mean()*3)/3 # to get mean value for codon
        else:
            sys.exit("-bkg can be 1 or 2  but is {}".format(args.bkg))
        
        if bkg <= thres_bkg:    #FILTER 2
            continue
        # -17,-47,-77,-107
        peak_stop = s.loc[-18 : -16].sum()    # codon coverage at STOP
        peak_30   = s.loc[-48 : -46].sum()    # 30 nt before 
        peak_60   = s.loc[-78 : -76].sum()    # 60 nt before
        peak_90   = s.loc[-108:-106].sum()    # 90 nt before
        # ratios - ivide codon coverage by avarage codon coverage  (rpm/rpm) - unitless
        qs_stop =  peak_stop/bkg 
        qs_30   =  (peak_30/bkg) 
        qs_60   =  (peak_60/bkg) 
        qs_90   =  (peak_90/bkg) 
        
        ##########
        # collect data
        c2+=1
        d[gtf.gene_id] = [bkg, qs_stop, qs_30, qs_60, qs_90]
        
print("{:8,} passed region threshold {}".format(c2,args.th1))
# pmake dataframe
index = ["bkg", "qs_stop", "qs_30", "qs_60", "qs_90"]
df = pd.DataFrame(d, index=index).T

# queuingScore-Weigthed
df['QS_weighted'] = df["qs_stop"] + df["qs_30"]*2 + df["qs_60"]*3 + df["qs_90"]*4
# queuingScore
df['QS']          = df["qs_stop"] + df["qs_30"] + df["qs_60"] + df["qs_90"]

if   weighted == "NO":
    df.drop(labels="QS_weighted", axis=1, inplace=True)
elif weighted == "YES":
    df.drop(labels="QS", axis=1, inplace=True)
elif weighted == "BOTH":
else:
    print("Warnings! parameter -w {} not recognised!".format(weighted))


i = 1   #  minimal queing score to keep in table
j = 50  #  strong queing score worth to check

mask1 = df['QS']>i
mask2 = df['QS']>j

df = df.loc[mask1,:]
print("There are  {} genes with 'queuingScore-W' bigger than {}".format(df.shape[0], i))
print("There are  {} most wavier genes {} times above bakground".format(df.loc[mask2,:].shape[0], j))
df.to_csv(outfile, sep='\t', header=True, index=True)
print("Output:\n   {}".format(outfile))
print("//")
hd5.close()
