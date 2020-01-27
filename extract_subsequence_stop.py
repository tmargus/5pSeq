#!/usr/bin/env python
#

__author__		= "Tonu Margus"
__copyright__	= "Copyright 2020"
__version__		= "0.1.5"
__email__		= "tonu.margus@gmail.com"
__status__		= "beta"

import sys
import argparse
import pandas as pd
import pysam


parser = argparse.ArgumentParser(description='Extracts sequence around STOP codon based on given gene names and threshold')
parser.add_argument('-annot',  type=str, help='GTF annotation file', default='0-References/genome.gtf.gz')
parser.add_argument('-genome',  type=str, help='Genome sequnece in FastA', default='0-References/Genome.fa')
parser.add_argument('-b',  type=int, help='Positions befoe stop in nt', default=30)
parser.add_argument('-a',  type=int, help='Positions after stop in nt', default=9)
parser.add_argument('-list',  type=str, help='filename with gene names in interest one per line', default=None)
parser.add_argument('-outpf',  type=str, help='Output format can be Fasta; Table;  plain', default='plain')
parser.add_argument('-translate',  type=str, help='Translate nt to aa sequence "Yes" or "No"', default='No')
parser.add_argument('-test',  type=str, help='print test message idf yes', default=None)

args = parser.parse_args()

sys.stderr.write("\n\
-annot   annotation GTF: {}\n\
-genome    genome FastA: {}\n\
-b    before stop codon: {}\n\
-a    after  stop codon: {}\n\
-translate    translate: {}\n\
-list   subsets by list: {}\n\
-outpf  subsets by list: {}\n".format(args.annot, args.genome, args.b, args.a,  args.translate,  args.list, args.outpf) )

usage = "./extract_subsequence_stop.py -annot 0-References/genome.gtf.gz  -genome  0-References/Genome.fa   -b 30 -a 9  > outfile.seq"

if args.test:
     sys.exit("\n  usage:\n\t{}\n".format(usage))
     
     
a = args.a         # nt after  stop codon
b = args.b         # nt before stop codon
annotation = args.annot
genomefile = args.genome
glist_file = args.list
translate  = args.translate.upper()
outputform = args.outpf.upper()

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
    return {info[0]: seq for info, seq in read_FASTA(filename)}

def get_part_from_gtf(annotation, reference=None, feature="CDS"):
    """ Returns a part from GTF annotation 
    annotation: "0-References/genome.gtf.gz"
    reference:  None # or chr in GTF file   'I' or 'XII' or 'chr1'
    feature:   'CDS' # or 'ORF' ar any valid feature in GTF    
    is file name to compressed and indexed GTF. vt tabix """
    
    tabixfile = pysam.TabixFile(annotation, parser=pysam.asGTF())
    return [gtf for gtf in tabixfile.fetch(reference=reference) if (gtf.feature == feature)]

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


genome = read_FASTA_dictionary(genomefile, SplitHeader=False)

stop_gtf_l = get_part_from_gtf(annotation, reference=None, feature="stop_codon") # gtf part for stop

# restrict by gene list if list is given
if glist_file is not None:
    try:
        # open file and read the content in a list
        with open(glist_file, 'r') as filehandle:  
            genes_l = [current_place.rstrip() for current_place in filehandle.readlines()]
            
        stop_gtf_l = [gtf for gtf in stop_gtf_l if gtf.gene_id in genes_l] 
    except:
        print('Could not open list file "{}"'.format(glist_file))
        print('Switching back to defaults (all genes)')

c = 0
l = [] 

for gtf in stop_gtf_l:
    c+=1            # count
    i = gtf.start   # start
    left_most, right_most = (i-b, i+a+3) if gtf.strand == '+' else (i-a, i+b+3) # define left & rightmost
    seq = genome[gtf.contig][left_most:right_most]      # get seq
    seq = seq if gtf.strand == '+' else revcompl(seq)   # reverese complement for - 
    seq = translate_DNA(seq) if translate == "YES" else seq # tranlsate if needed
        
    l.append(seq)
    # stdout
    if outputform == "FASTA":
        print(">{}\n{}".format(gtf.gene_id, seq))
    elif outputform == "TABLE":
        print("{}\t{}".format(gtf.gene_id, seq))
    else:
        print(seq)

sys.stderr.write("\nsequences: {}\n".format(c) )

