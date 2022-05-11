# Contig stitcher

A script in Python 3 that stitches contigs from sequencing together to create pseudo-references files for downstream sequencing analysis. 

## Description

Contig stitcher is developed in Python 3, and had been confirmed to run on a Mac OSX environment (not installation needed!). 

Contig stitcher takes a supplied fasta (.fa/.fna) file containing multiple sequences, and concatenates them together using spacers of 100 N's to generate a new output file. Users can specify the maximum length of these 'pseudocontigs', in which case multiple pseudocontigs will be generated in the output file. Contig stitcher also generates a .bed file showing the location of each original sequence within the pseudocontig(s). 

## Usage

```
usage: contig-stitcher.py [-h] [-l LENGTH] [-v] [-n] reference

positional arguments:
  reference             The main fa/ fna file to be concatenated

options:
  -h, --help            show this help message and exit
  -l LENGTH, --length LENGTH
                        The max length of each pseudocontig. Default is 1 Mbp
  -v, --verbose         Print additional commentary
  -n, --numbers         Use only numbers for pseudocontig labels
```
