#!/usr/bin/env python
# 
# Script for 5' assignment of 5'P-Seq data
# input is BAM file  must contain NH tag
# reads with the tag NH:i:1   only included
# output 1: raw counts in *_iv.h5          - single indexed
# output 2: normalised RPM in _idx_iv.h5   - double indexed
#
__author__		= "Tonu Margus"
__copyright__	= "Copyright 2020"
__version__		= "0.1.2"
__email__		= "tonu.margus@gmail.com"
__status__		= "beta"

import re
import sys
import pysam
import argparse
import pandas as pd
from collections import Counter

parser = argparse.ArgumentParser(description="five prime assinment of 5'P-Seq data" )
parser.add_argument('-i',  type=str, help='aligned_sorted.bam')
args = parser.parse_args()

sys.stderr.write("\n\
-i  input           : {}\n\n".format(args.i))

usage = 'python fivePassignment.py -i aligned_sorted.bam"'

if args.i==None:
     sys.exit("\n  usage:\n\t{}\n".format(usage))

raw_out = False # bool 

# output file name from infilename
f = re.sub(r'_sorted.bam', '', re.sub(r'.*\/', '', args.i))
outf_raw_hdf = "{}_raw_iv.h5".format(f)
outf_rpm_hdf = "{}_rpm_iv.h5".format(f)
outf_idx_hdf = "{}_idx_iv.h5".format(f)

def yeastChr():
    # Ordered yeast Chr list short names from ensembl
    return ['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI','Mito']

def update_df(df, Chr, strand):
    df.fillna(0, inplace=True)
    columns = list(df.columns)
    columns = ["Chr", "Position", "Strand"] + columns
    df["Chr"] = Chr
    df["Strand"] = strand
    df["Position"] = df.index

    return df[columns]

def restructurate_hd5(infile, outfile, close_outfile=True):
    """ infile.h5 keys - "/For_raw", "/Rev_raw", ...
    outfile.h2  keys - "/For_raw/I", "/For_raw/II", ... etc
    "Position" is set to index

    :param infile:
    :param outfile:
    :return: reindexed 2 level hdf
    """
    # open inp_HDF
    inp__h5 = pd.HDFStore(infile, "r")
    outp_h5 = pd.HDFStore(outfile, complevel=5, complib="zlib", mode="w")
    # open out_HDF

    # for each I level table For_raw, Rev_raw, ...
    for key in inp__h5.keys():
        # for each chromosome
        df = inp__h5[key]
        for Ch in df['Chr'].unique():
            df_ch = df[df.Chr == Ch].copy()

            # set Position to index
            df_ch.set_index('Position', inplace=True)

            # save df under II level key what is now chromosome
            h5_key = key + "/" + Ch
            outp_h5.put(h5_key, df_ch)

    inp__h5.close()

    if close_outfile == True:
        outp_h5.close()
    else:
        return outp_h5

    outp_h5.close()


# empty dataframe for collecting data
df_for_sum = pd.DataFrame()
df_rev_sum = pd.DataFrame()

bamfile    = pysam.AlignmentFile(args.i, "rb")  # open BAM file
sys.stderr.write("reads mapped once:\n")
for ref in yeastChr():
    # assumes NH:i:  tag in BAM file -> counts reads mapped once
    # list of forward strand reads 5' coordinates
    For_l = [read.reference_start for read in bamfile.fetch(ref) if (not read.is_reverse & read.get_tag("NH")==1)]
    For_s = pd.Series(Counter(For_l))
    For_s.name = "counts"
    # list of reverse strand reads 5' coordinates
    Rev_l = [read.reference_end - 1 for read in bamfile.fetch(ref) if (read.is_reverse & read.get_tag("NH")==1)]
    Rev_s = pd.Series(Counter(Rev_l))
    Rev_s.name = "counts"
    report = "{:4} {:>9,}\n".format(ref, len(For_l)+len(Rev_l))
    sys.stderr.write(report)
    df_for = update_df(pd.DataFrame(For_s), Chr=ref, strand="+")
    df_rev = update_df(pd.DataFrame(Rev_s), Chr=ref, strand="-")

    df_for_sum = pd.concat([df_for_sum, df_for], ignore_index=True)  # collect summary table
    df_rev_sum = pd.concat([df_rev_sum, df_rev], ignore_index=True)  # collect summary table

# raw counts
if raw_out:
    store = pd.HDFStore(outf_raw_hdf, complevel=5, complib="zlib", mode="w")
    store.put("For_raw", df_for_sum, format="table", data_columns=True)
    store.put("Rev_raw", df_rev_sum, format="table", data_columns=True)
    store.close()
    sys.stderr.write("Raw counts\t    : {}\n".format(outf_raw_hdf))
else:
    pass
    
# raw -> rpm
l = [0 for read in bamfile.fetch() if read.get_tag("NH") == 1]  # reads mapped once
normFactor = len(l) / 10 ** 6  # normalisation factor
report = "normalization factor: {}\nreads mapped once   : {:,}\n".format(normFactor, len(l))

columns = ["Chr","Position","Strand","rpm"]
df_for_sum["rpm"] = df_for_sum["counts"] / normFactor
df_rev_sum["rpm"] = df_rev_sum["counts"] / normFactor

store = pd.HDFStore(outf_rpm_hdf, complevel=5, complib="zlib", mode="w")
store.put("For_rpm", df_for_sum[columns], format="table", data_columns=True)
store.put("Rev_rpm", df_rev_sum[columns], format="table", data_columns=True)
sys.stderr.write("Normalised rpm\t    : {}\n".format(outf_rpm_hdf))
store.close()
bamfile.close()

# restructureate
restructurate_hd5(outf_rpm_hdf, outf_idx_hdf, close_outfile=True)
sys.stderr.write("Restructurate \t    : {}\n".format(outf_idx_hdf))
sys.stderr.write(report)