# hpSource
## Purpose
A set of tools to find, analyze & visualize hairpin RNAs (hpRNAs)
## hpSource.py
Finds hairpin RNAs in a genome and generates annotation file in bed format.

Written in python3.  

Dependencies:

[biopython](https://github.com/biopython/biopython.github.io/)

[einverted](http://emboss.bioinformatics.nl/cgi-bin/emboss/help/einverted)

[revseq](http://www.bioinformatics.nl/cgi-bin/emboss/help/revseq)

[bedtools](http://bedtools.readthedocs.io/en/latest/)

RNAfold from the [ViennaRNA](https://www.tbi.univie.ac.at/RNA/#download) package

Basic usage is:
```bash
hpSource.py -g input.fas -s sRNA.fastq
```
hpSource takes two mandatory arguments:

	-g (input genome fasta file)

	-s (input small RNA fastq file)

hpSource takes five optional arguments:

	-p (gap penalty for einverted, default=12)

	-t (minimum score threshold for einverted, default=50)

	-m (match score for einverted, default=3)

	-a (mismatch score for einverted, default=-4)

	-f (threshold value for fraction of coverage contributed by dominant strand)
