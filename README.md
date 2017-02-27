# hpSource
## Purpose
A set of tools to find, analyze & visualize hairpin RNAs (hpRNAs)
## hpSource.py
Finds hairpin RNAs in a genome and generates gff annotation.

Written in python3.  

Dependencies:

[einverted](http://emboss.bioinformatics.nl/cgi-bin/emboss/help/einverted)

Basic usage is:
```bash
hpSource.py -i input.fas
```
hpSource takes one mandatory argument:

	-i (input fasta file)

hpSource takes four optional arguments:

	-g (gap penalty for einverted, default=12)

	-t (minimum score threshold for einverted, default=50)

	-m (match score for einverted, default=3)

	-s (mismatch score for einverted, default=-4)

