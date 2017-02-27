#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os
import re

# argument parsing
parser=argparse.ArgumentParser(description="Read arguments")
parser.add_argument("-g","--inputgenome",type=str,help="Input genome file")
parser.add_argument("-s","--inputsrna",type=str,help="Input sRNA read file")
parser.add_argument("-p","--gap",type=int,help="Gap penalty (default=12)")
parser.add_argument("-t","--threshold",type=int,help="Minimum score threshold (default=50)")
parser.add_argument("-m","--match",type=int,help="Match score (default=3)")
parser.add_argument("-a","--mismatch",type=int,help="Mismatch score (default=-4)")
args=parser.parse_args()
# default values
EinvGap=12
EinvThreshold=50
EinvMatch=3
EinvMismatch=-4
# parse input genome file
if args.inputgenome is not None:
	InputGenome = args.inputgenome
	print("Input genome = "+InputGenome)
else:
	print("ERROR: input genome file (-g) not specified")
	sys.exit(0)
# parse input sRNA file
if args.inputsrna is not None:
	InputsRNA = args.inputsrna
	print("Input sRNA file = "+InputsRNA)
else:
	print("ERROR: input sRNA file (-s) not specified")
	sys.exit(0)
# parse gap penalty
if args.gap is None:
	print("Using default gap penalty ("+str(EinvGap)+")")
else:
	if args.gap >= 0:
		EinvGap=int(args.gap)
		print("Gap penalty = "+str(EinvGap))
	else:
		print("ERROR: gap penalty (-g) must be an integer of 0 or more")
		sys.exit(0)
# parse minimum score threshold
if args.threshold is None:
	print("Using default minimum score threshold ("+str(EinvThreshold)+")")
else:
	if args.threshold >= 0:
		EinvThreshold=int(args.threshold)
		print("Minimum score threshold = "+str(EinvThreshold))
	else:
		print("ERROR: minimum score threshold (-t) must be an integer of 0 or more")
		sys.exit(0)
# parse match score
if args.match is None:
	print("Using default match score ("+str(EinvMatch)+")")
else:
	if args.match >= 0:
		EinvMatch=int(args.match)
		print("Match score = "+str(EinvMatch))
	else:
		print("ERROR: match score (-m) must be an integer of 0 or more")
		sys.exit(0)
# parse mismatch score
if args.mismatch is None:
	print("Using default mismatch score ("+str(EinvMismatch)+")")
else:
	if args.mismatch < 0:
		EinvMismatch=int(args.mismatch)
		print("Mismatch score = "+str(EinvMismatch))
	else:
		print("ERROR: mismatch score (-s) must be an integer below 0")
		sys.exit(0)

for i in [".fa",".fas"]:
	if InputGenome.endswith(i):
		os.rename(InputGenome,InputGenome.replace(i,".fasta"))
		InputGenome=InputGenome.replace(i,".fasta")

# function to run revseq (to complement genome, to find antisense hpRNA)
def revseqRun(genome):
	print("Generating complement of "+genome)
	cmd="revseq -sequence "+genome+" -outseq "+genome.replace(".fasta","_RC.fasta")
	subprocess.call(cmd,shell=True)
	print("Reverse complement of sequence "+genome+" written to "+genome.replace(".fasta","_RC.fasta"))
	return(genome.replace(".fasta","_RC.fasta"))

# function to run einverted
def einvertedRun(genome):
	print("Finding secondary structure in "+genome) 
	cmd="einverted -sequence "+genome+" -gap "+str(EinvGap)+" -threshold "+str(EinvThreshold)+" -match "+str(EinvMatch)+" -mismatch "+str(EinvMismatch)+" -outfile "+genome.replace(".fasta",".out")+" -outseq "+genome.replace(".fasta","_einverted.fasta")
	subprocess.call(cmd,shell=True)
	print("Output of einverted written to "+genome.replace(".fasta",".out")+" & "+genome.replace(".fasta","_einverted.fasta"))
	return(genome.replace(".fasta",".out"))

# function to convert einverted output to gff
def einvertedParse(results):
	Locations=[]
	Names=[]
	Scores=[]
	LeftStarts=[]
	LeftEnds=[]
	RightStarts=[]
	RightEnds=[]
	print("Parsing einverted output file: "+results)
	chrom=''
	LeftArm=True
	hpCount=1
	for line in open(results,"r"):
		if line.startswith("         "):
			LeftArm=False
		elif re.search("\:",line):
			score=int(line.split(" Score ")[1].split(":")[0])
			Scores.append(score)
			if chrom==line.split(":")[0]:
				hpCount+=1
				Locations.append(chrom)
				Names.append(chrom+"_hpRNA"+str(hpCount))
			else:
				chrom=line.split(":")[0]
				hpCount=1
				Locations.append(chrom)
				Names.append(chrom+"_hpRNA"+str(hpCount))
			LeftArm=True
		elif line.startswith("   "):
			print(line)
			temp=line.strip("   ").split(" ")
			if LeftArm==True:
				start=int(temp[0])
				LeftStarts.append(start)
				end=int(temp[2])
				LeftEnds.append(end)
			elif LeftArm==False:
				start=int(temp[2])
				RightStarts.append(start)
				end=int(temp[0])
				RightEnds.append(end)
	print(str(len(Locations))+" locations found")
	print(str(len(Names))+" names found")
	print(str(len(Scores))+" scores found")
	print(str(len(LeftStarts))+" left starts found")
	print(str(len(LeftEnds))+" left ends found")
	print(str(len(RightStarts))+" right starts found")
	print(str(len(RightEnds))+" right ends found")
	BedOutput=""
	for i in range(len(Locations)):
		BedOutput+=Locations[i]+"\t"+str(LeftStarts[i])+"\t"+str(LeftEnds[i])+"\t"+Names[i]+"_left\t"+str(Scores[i])+"\n"+Locations[i]+"\t"+str(RightStarts[i])+"\t"+str(RightEnds[i])+"\t"+Names[i]+"_right\t"+str(Scores[i])+"\n"
	output = open(results.replace(".out",".bed"),"wt")
	output.write(BedOutput)
	output.close()
	print("Bed file of hpRNAs written to "+results.replace(".out",".bed"))
	return()

# function to map sRNAs to a genome
def sRNAmap(srna,genome):
	if srna.endswith(".gz"):
		print("Input file "+srna+" is compressed - uncompressing...")
		cmd="gunzip "+srna
		subprocess.call(cmd,shell=True)
		print("Input file "+srna+" uncompressed")
		srna=srna.replace(".gz","")
	# make bowtie2 database for input genome
	print("Making bowtie2 database for genome "+genome)
	cmd="bowtie2-build "+genome+" ./"+genome.replace(".fasta","")
	subprocess.call(cmd,shell=True)
	# map reads to forward strand
	print("Mapping "+srna+" to forward strand of "+genome)
	cmd="bowtie2 --fast --norc -q -x "+genome.replace(".fasta","")+" -U "+srna+" -S MappedF.sam"
	subprocess.call(cmd,shell=True)
	print("Reads mapped")
	print("Converting SAM to BAM")
	cmd="samtools view -bS -u MappedF.sam > MappedF.bam"
	subprocess.call(cmd,shell=True)
	print("SAM converted to BAM")
	os.remove("MappedF.sam")
	# map reads to reverse strand
	print("Mapping "+srna+" to reverse strand of "+genome)
	cmd="bowtie2 --fast --nofw -q -x "+genome.replace(".fasta","")+" -U "+srna+" -S MappedR.sam"
	subprocess.call(cmd,shell=True)
	print("Reads mapped")
	print("Converting SAM to BAM")
	cmd="samtools view -bS -u MappedR.sam > MappedR.bam"
	subprocess.call(cmd,shell=True)
	print("SAM converted to BAM")
	os.remove("MappedR.sam")
	return()

# master function to: run einverted
def hpRNAfind(input):
	inputRC=revseqRun(genome=input)
	einvertedResults=einvertedRun(genome=input)
	einvertedParse(results=einvertedResults)
	einvertedResultsRC=einvertedRun(genome=inputRC)
	einvertedParse(results=einvertedResultsRC)	
	return("hpRNA analysis complete for "+input)

hpRNAfind(input=InputGenome)
#sRNAmap(srna=InputsRNA,genome=InputGenome)
