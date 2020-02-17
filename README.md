# 5pSeq
Scripts for 5' Seq data analysis in Python 3 (ver 3.6)are split in two groups based where reads were mapped:  a) to genome, b) to cDNA.  
Scripts are in development stage so don't use for publications.
In addition to basic Python libraries it assumes `pysam`, `argparse`, `pandas`, `collections`, `glob` to be installed.   

## Reads mapped to GENOME
```bash
# data files in the folder    0-References/
Genome.fa # is link to Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fasta
Genome.gtf # is link to       Saccharomyces_cerevisiae.R64-1-1.95.gtf
#
# reformated and indexed annotation for tabix (part of pysam)
genome.gtf.gz
genome.gtf.gz.tbi
```

### Mapping 5PSeq data
`fivePSeqMap.py`  
maps reads 5' ends. Output is in compressed HDF5 format. HDF file have keys for each strand chromosome; for example: `"For_rpm/I"`   and    `"Rev_rpm/I"`  for the first chromosome.

```bash
# Run in batch
for f in 5-Aligned/*.bam; do ./fivePSeqMap.py -i $f; done 
```

### Metagene plots  before STOP
1. Get coverage around stop  
	`hdf_2_metagene_tables_stop_5PSeq.py`
Creates metagene tables for stop codon. Input is output of previous script (`fivePSeqMap.py` )   files ending with `*_idx_iv.h5`.  Output is table with relative coordinates from -120  (default value in nt inside gene) up to 120. 0 position corresponds to 5' position of stop codon.
```bash
# command line parameters 
./hdf_2_metagene_tables_stop_5PSeq.py 

-i      input *.h5:      None
-prefix part of oufile : None
-annot   annotation GTF: 0-References/genome.gtf.gz
-genome    genome FastA: 0-References/Genome.fa
-col            columns: rpm
-th    region threshold: 15
-span   included region: 120
-norm  to gene coverage: Yes
-subsets  by stop codon: None
-glist  subsets by list: None


  usage:
	./hdf_2_metagene_tables_1_dev.py  -i WTS1_5-End_21-33_idx_assign_rpm.h5   -prefix WTS1_stop_metagene   -norm  Yes
```

```bash
# batch run
for f in 6-AssignRaw/*.h5; do fo=${f##*/}; ./hdf_2_metagene_tables_stop_5PSeq.py -i $f -prefix ${fo/_idx_iv.h5} -col rpm  -th 18 -span 180 -subsets NO ; done | tee hdf_2_metagene.log
```

2. Merge coverage of different samples
`merge_metagene_tbl_5PSeq.py`

```bash
./merge_metagene_tbl_5PSeq.py -files "8-MetagTbl/*sum*" -o meta-Stop_th18-Span180.csv

```
3. Plot metagene  
	 `Plot-STOP-MetaG-5PSeq.ipynb`
	 
### Computing queuing score
`compute_queuing_5PSeq.py` 
Computes weighted queuing for each gene with a coverage above threshold.
```bash
./compute_queuing_5PSeq.py -h

Computes ribosomes queing score at stop

  -h, --help    show this help message and exit
  -i I          input table of gene coverage in *.hd5 format
  -o O          output file name *.csv
  -annot ANNOT  GTF annotation file
  -th1 TH1      Summary gene coverage 150 nt before stop - 10(rpm) default
  -th2 TH2      Background coverage - codon mean from -115 up to Span
  -span SPAN    Positions before - stop recommended 150 or bigger
  -bkg BKG      region for bakground con be 1 or 2: 1 - before and 2 - between
                peaks
  -col COL      column for values: "sum"; "rpm"; "counts"
```

```bash
for f in 6-AssignRaw/*.h5; do fo=${f##*/}; ./compute_queuing_5PSeq.py -i $f -th1 15 -th2 0.15 -span 180 -col rpm -o ${fo/_idx_iv.h5/_queuingScores.csv}; done
```

### Extract sequence around stop codon
`extract_subsequence_stop.py`

```bash
./extract_subsequence_stop.py -test Yes

-annot   annotation GTF: 0-References/genome.gtf.gz
-genome    genome FastA: 0-References/Genome.fa
-b    before stop codon: 30
-a    after  stop codon: 9
-translate    translate: No
-list   subsets by list: None
-outpf  subsets by list: plain

  usage:
	./extract_subsequence_stop.py -annot 0-References/genome.gtf.gz  -genome  0-References/Genome.fa   -b 30 -a 9  > outfile.seq


# last amino acid for all genes
./extract_subsequence_stop.py -b 3 -a -3 -translate Yes -outpf Table  > last_aa_tbl.txt

# last coding codon
./extract_subsequence_stop.py -b 3 -a -3  -outpf Table  > last_coding-codon_tbl.txt

# stop codon
./extract_subsequence_stop.py -b 0 -a 0  -outpf Table  > stop-codon_tbl.txt

# stop codon extended by one towards  3 prime
./extract_subsequence_stop.py -b 0 -a 1  -outpf Table  > stop-codon-extended-3pr_tbl.txt

# stop codon extended by one towards  5 prime
./extract_subsequence_stop.py -b 1 -a 0  -outpf Table  > stop-codon-extended-5pr_tbl.txt
```


## Reads mapped to cDNA

```bash
# data file
Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa
```

Main thing to consider is that cDNA sequences  in most cases are starting with AUG (start codon)  and ending with stop codon. Therefore nothing before start codon is aligned and region before stop is also poorly covered. Reads are usually 60 nt long and if half of read don't align (because there is no sequence towards 3' end in this cDNA sequence) read is not aligned. This setup complements to genome. 
 
### Map 5' positions cDNA  
`fivePSeqMap2cdnawithOffset.py`
Input is:  
1.  cDNA aligned reads in sorted indexed BAM and 
2. fastA file of cDNAs
Output: offset (17 A-Site) corrected coverage for each cDNA
1. pickled dictionary with cDNA IDs as keys and pandas table as value
Tables index corresponds to codon positions. Columns are "codon"; "amino acid" and "counts" 

```bash
usage: fivePSeqMap2cdnawithOffset.py [-h] [-i I] [-f F] [-of OF] [-o O]

five prime assinment (offset 17) of cDNA aligned 5'P-Seq reads

optional arguments:
  -h, --help  show this help message and exit
  -i I        aligned_to_cDNA_sorted.bam
  -f F        Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa
  -of OF      offset for correcting 5prim mapping; default 17 (A-Site)
  -o O        output pickle file
```

```bash
 #P-site offset 15
fivePSeqMap2cdnawithOffset.py -i 5-cdna-Aligned/5PSeq/KO20C_V_5PSeq_dedup.bam -of 15 -o KO20C_V_5PSeq_cdna_off15.pkl
```


### Compute QS for fixed codon in the middle of a gene

```bash
fivePSeqQSfromcDNA4PSiteCodon.py  -codon AAA  -i KO20C_V_5PSeq_cdna_off15.pkl -o KO20C_V_QS-coding-region_AAA_P-site.csv

-i     input pickled dict : KO20C_V_5PSeq_cdna_off15.pkl
-codon     P-Site codon   : AAA
-per periodicity in codons: 10
-o     output_table.csv   : KO20C_V_QS-coding-region_AAA_P-site.csv

4,905 cDNAs from input file
QS for P-Site AAA  A-Site NNN pairs
10,000 ...
20,000 ...
30,000 ...
40,000 ...
50,000 ...
60,000 ...
70,000 ...
80,000 ...
83,534 total QS for codon pairs AAA[NNN]!
KO20C_V_QS-coding-region_AAA_P-site_v2.csv

```
