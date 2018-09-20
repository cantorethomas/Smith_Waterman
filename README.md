# Smith_Waterman

## What does it do
This is an implementation of Smith Waterman algorithm for local alignments. 
In this version, it is possible to print out the first *n-highest scores* alignments for each couple of sequences in input.\
User can also set specific thresholds for minimum length and score of the alignments.\
Match, mismatch and gap score can be customized: default are respectively +3, -3 and -2. By the way, default one can be easily
modified at the beginning of the python script.

## Options list 

```
--seq1          The first sequence to align
--seq2          The second sequence to align
-i, --input     input file 
-m, --match     Match score
-s, --mismatch  Mismatch score 
-g, --gap       Gap score
--minscore      Minimum score
--minlength     Minimum length 
--numresult     Number of alignments to return 
```

## Input format
It is possible to give sequences in input directly from the command line using options `--seq1` and `--seq2`, or by giving a file
in input, each line consisting of two sequences in a row, tab separated.

```
# three couple of sequences to align 
ATCGATGCGGAGCTTAGGCTGA  AGCGCGATGGATGAGAGGCTGAT
AGCTGAGCTGGAGGGGAGTCGA  AGGGCGACGCGATGACGCACCAG
AGCGTGGATGGAGGAGTGTTTT  AGTGAGTGTGACGCGCGATTTAC
```
## For a quick example 
"Example" folder contains the simple `input_gen.py` script that generates a random input file for a quick test of the main program. The input 
file can be generated running the script from bash choosing the length of each sequence (`--nseq`) and the total number of couple of sequences
to be aligned (`--len`).

