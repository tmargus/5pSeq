# 5pSeq
Scripts for 5' Seq data analysis in Python 3 (ver 3.6).
Scripts are in development stage so don't use for publications.

In addition to basic Python libraries it assumes `pysam`, `argparse`, `pandas`, `collections`, `glob` to be installed. 

## Short description of scripts
### Mapping 5PSeq data
`fivePSeqMap.py`   
maps reads 5' ends. Output is in compressed HDF5 format. 
```bash
# Run in batch
for f in 5-Aligned/*.bam; do ./fivePSeqMap.py -i $f; done 

```
### Metagene plots  before STOP
1. Get coverage around stop  
`hdf_2_metagene_tables_stop_5PSeq.py`
Creates metagene tables for stop codon. Input is output of previous script (`fivePSeqMap.py` )   files ending with `*_idx_iv.h5`.  Output is table with relative coordinates from -120  (default value in nt inside gene) up to 120. 0 position corresponds to 5' position of stop codon.  
```bash
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

optional arguments:
  -h, --help        show this help message and exit
  -i I              input table of gene coverage in input.hd5 
  -o O              output file name output_tbl.csv
  -annot ANNOT      GTF annotation file
  -th1 TH1          Summary gene coverage 180 nt before stop - 10(rpm) default
  -th2 TH2          Background coverage - codon mean from -115 up to Span
  -span SPAN        Positions before - stop recommended 180 or bigger
  -col COL          column for values: "sum"; "rpm"(default); "counts"
  -subsets SUBSETS  Split to subsets based Stop codon
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
