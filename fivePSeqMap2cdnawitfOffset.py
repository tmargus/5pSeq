#!/usr/bin/env python
# Python 3
#
# Script for 5' assignment of 5'P-Seq data
# input is aligned reads to cDNA 
# for example: Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa 
# BAM file is sorted and indexed with reads mapped once, i. e. the tag NH:i:1
# Output:
# is dictionary of pandas tables: 
#  with the cDNA ID as  a key 
#  each table come with the index and three columns: 
#  index: position of codon
#  codon: 
#  amino acid: single code
#  coverage: 5' end of reads mapped  and summed up per codon
#  reads are already corrected by offset. Ofset 17 - A-Site; OF 15 P-Site

__author__		= "Tonu Margus"
__copyright__	= "Copyright 2020"
__version__		= "0.1"
__email__		= "tonu.margus@gmail.com"
__status__		= "beta"

import re
import sys
import pysam
import argparse
import pandas as pd
import pickle
from collections import Counter

parser = argparse.ArgumentParser(description="five prime assinment (offset 17) of cDNA aligned 5'P-Seq reads" )
parser.add_argument('-i',  type=str, help='aligned_to_cDNA_sorted.bam')
parser.add_argument('-f',  type=str, help='Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa', default='0-References/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa')
parser.add_argument('-o',  type=str, help='output pickle file', default='output_cdna.pkl')
args = parser.parse_args()

## set some filters
min_counts   = 10  # minimal number of reads mapped to gene
min_aver_cov = 0.2 # minimal average coverage  all_counts/gene_length_in_codons
 ##
sys.stderr.write("\n\
-i  input sorted BAM : {}\n\
-f  cDNA in fasta    : {}\n\
-o  output pkl       : {}\n\n".format(args.i, args.f, args.o))

usage = 'python fivePSeqMap2cdnawitfOffset.py  -i aligned_sorted.bam -o  output_cdna.pkl"'

if args.i==None:
     sys.exit("\n  usage:\n\t{}\n".format(usage))
###################################################################
# Params and I/O
offset = 17 # 5' offset: 17 is A-site correction
pickle_f = args.o
bam = args.i
bamfile = pysam.AlignmentFile(bam, "rb")  # open BAM file
# for data collection
gene_codon_d = {}
# counters
c = cc = 0
low_coverage = 0
low_average_coverage = 0
# 
print("Mapping reads to cDNA!\nOnly Forward reads are counted!")
print("offset = {}".format(offset)) 
###################################################################
# Functions
def read_FASTA(filename, splitstr='|', SplitHeader=True):
    """ Reads FastA file and returns a list of tuples, where first 
    part is a list of header elements and second seq as a string 

    read_FASTA('seqfile.fa', SplitHeader=True) 
        [(['gi', '1114567', 'gb', 'NC_00245'],
        'ATATAGGCGCTTGGTGCGCGGCGGGCGCGGCTAGCAGCACCTTTAGTAGCTTTCATCAT'),
        (['gi', '2224567', 'gb', 'NC_22245'],
        'ATTTTTGGGGGGCGCGGCTAGCAGCACCTTTAGTAGCTTTCAAAAAATTTTCAT')]
    
    info is:
        >gi|1114567|gb|NC_00245

    "Bioinformatics Programming Using Python by Mitchell L Model"
    """
    with open(filename) as file:

        if SplitHeader:
            return [(part[0].split(splitstr),
                part[2].replace('\n', ''))
                for part in 
                [entry.partition('\n')
                    for entry in file.read().split('>')[1:]]]
        else: 
            return [(part[0],
                part[2].replace('\n', ''))
                for part in 
                [entry.partition('\n')
                    for entry in file.read().split('>')[1:]]]

def read_FASTA_dictionary(filename, splitstr='|', SplitHeader=False):
    """ Creates dictionary from FastA file, where key is gi-number and value is seq 

    make_indexed_sequence_dictionary('seqfile.fa')
        {'1114567': 'ATATAGGCGCTTGGTGCGCGGCGGGCGCGGCTAGCAGCACCTTTAGTAGCTTTCATCAT',
         '2224567': 'ATTTTTGGGGGGCGCGGCTAGCAGCACCTTTAGTAGCTTTCAAAAAATTTTCAT' }
    
    read_FASTA by default splits header '|' assuming NCBI entry but 
    here read_FASTA do not split header  (SplitHeader=False), i.e. key is the whole name.
    
    "Bioinformatics Programming Using Python by Mitchell L Model"
    """
    return {info[0]: seq for info, seq in read_FASTA(filename, splitstr=splitstr, SplitHeader=SplitHeader)}

def complement(seq):
    basecomplement = {'A':'T','C':'G','G':'C','T':'A','Y':'R','R':'Y','M':'K','K':'M','W':'W','V':'B','B':'V','H':'D','D':'H','N':'N',
                      'a':'t','c':'g','g':'c','t':'a','y':'r','r':'y','m':'k','k':'m','w':'w','v':'b','b':'v','h':'d','d':'h','n':'n' }
    return ''.join([basecomplement[base] for base in seq])

def revcompl(seq):
    return complement(seq[::-1])

# Translates 
def translate_DNA_codon(codon):
    ''' Translates DNA triplet to amino acid (20) by using single 
    lettr code where underline "_" corresponds to STOP codon.
    '''
    return DNA_codon_table[codon] 

def aa_generator_DNA(dnaseq):
    """Return a generator object that produces an amino acid by translating 
    the next three characters of dnaseq each time next is called on it"""
    return (translate_DNA_codon(dnaseq[n:n+3])
            for n in range(0, len(dnaseq), 3))

def translate_DNA(dnaseq):
    """Translate dnaseq into amino acid symbols"""

    gen = aa_generator_DNA(dnaseq)
    seq = ''
    aa = next(gen, None)
    while aa:
        seq += aa
        aa = next(gen, None)
    return seq

DNA_codon_table = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}

def update_df(df, cdna, strand):
    df.fillna(0, inplace=True)
    #df["sum"] = df.sum(axis=1)
    columns = list(df.columns)
    columns = ["cdna", "Position", "Strand"] + columns
    df["cdna"] = cdna
    df["Strand"] = strand
    df["Position"] = df.index

    return df[columns]

def series_2_df(data, offset=0):
    ''' data = [s1,s2]
    s1: series from cDNA  seq_s = pd.Series(list(seq),name='seq')
    s2: series of positions with some counts  For_s   = pd.Series(Counter(For_l)); For_s.name = "counts"
    There can be more than 2 series. Series must come with names
    '''
    df = pd.concat(data, axis=1)
    df.fillna(0, inplace=True)
    df['counts'] = df['counts'].shift(periods=offset, fill_value=0)       
    return df

def convert_nt_2_codon_based(df):
    '''Converts nucleotide based df to codon based df
    index will corresponf to codon #
    input comes from def series_2_df(data):
    output table contains columns: 
    'codon','animo acid','counts'
    '''
    d = {}
    seq_max = df.shape[0] # jääk %3 == 0
    
    for n,i in enumerate(list(range(0,seq_max,3))):
        codon = df_f1.loc[i:i+2,'seq'].sum()
        aa    = translate_DNA_codon(codon)
        count = df_f1.loc[i:i+2,'counts'].sum()
        d[n]  = (codon,aa,int(count))

    return pd.DataFrame(d, index=['codon','animo acid','counts']).T
    

###################################################################
# main program
fasta_file = args.f
fasta = read_FASTA_dictionary(fasta_file, splitstr=' ', SplitHeader=True)
keys = list(fasta.keys())
print("{:,} sequences  to be processed".format(len(keys)))

for ref in keys: #[k for k in keys if k  in included]:    
    if len(ref)<4:
        print("WARNINGS! Skipping unusual cDNA ID {}".format(ref))
        continue 
    seq = fasta[ref]
    if len(seq)%3 >0:
        print("WARNINGS! Skipping {}! cDNA is out of frame! length of cDNA is {:,} nt.".format(ref, len(seq)))
        continue
    if not seq:
        print("WARNINGS! Skipping {}! No sequence retrieved!!!".format(ref, len(seq)))
        continue
        
    c+=1
    if c == 200:
        cc += c
        c=0
        print("{:6,} ...".format(cc))        
    c1  = c2 = 0
    
    # list of forward strand reads 5' coordinates
    For_l = [read.reference_start for read in bamfile.fetch(ref) if not read.is_reverse]
    if not For_l:
        continue
    else:
        # get sequence
        
        seq_s   = pd.Series(list(seq),name='seq')
        For_s   = pd.Series(Counter(For_l), name='counts')
        if For_s.sum() < min_counts:     # sckip genes with less than 10 counts 
            low_coverage+=1
            continue
        if For_s.sum()/(len(seq)/3) < min_aver_cov:     # count genes with average low coverage 
            low_average_coverage+=1
            continue
            
        data    = [seq_s,For_s]
        df_f1   = series_2_df(data, offset=offset) # series to dataframe
        df_f1_codon = convert_nt_2_codon_based(df_f1) # convert to codon based table
        df_f1_codon['counts'] = df_f1_codon['counts'].apply(pd.to_numeric) # make counts numeric

        # add to dictionary
        gene_codon_d[ref] = df_f1_codon


print("{:,}\t # skipped low coverage ({}) ".format(low_coverage, min_counts))
print("{:,}\t # skipped low awerage coverage ({})".format(low_average_coverage, min_aver_cov))
print("{:,}  sequences included".format(len(gene_codon_d)))

with open(pickle_f , 'wb') as handle:
    pickle.dump(gene_codon_d, handle, protocol=pickle.HIGHEST_PROTOCOL)

print(pickle_f)
