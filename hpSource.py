#!/usr/bin/env/ python3

import argparse
import subprocess
import sys
import os

# argument parsing
parser=argparse.ArgumentParser(description="Read arguments")
parser.add_argument("-i","--input",type=str,help="Input genome file")
parser.add_argument("-g","--gap",type=int,help="Gap penalty (default=12)")
parser.add_argument("-t","--threshold",type=int,help="Minimum score threshold (default=50)")
parser.add_argument("-m","--match",type=int,help="Match score (default=3)")
parser.add_argument("-s","--mismatch",type=int,help="Mismatch score (default=-4)")
args=parser.parse_args()
# default values
EinvGap=12
EinvThreshold=50
EinvMatch=3
EinvMismatch=-4
# parse input genome file
if args.input is not None:
	InputGenome = args.input
	print("Input genome = "+InputGenome)
else:
	print("ERROR: input genome file (-i) not specified")
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
# make output name stem
#OutputStem=""
#for i in [".fa",".fas",".fasta"]:
#	if InputGenome.endswith(i):
#		OutputStem=InputGenome.replace(i,"")
#		print("Output stem name set to "+OutputStem)
#if OutputStem=="":
#	print("ERROR: output stem could not be created. If input fasta file does not end with .fa or .fas or .fasta, rename it to one of these and try again")
#	sys.exit(0)

# if input genome file ends in .fa or .fas, rename it to end in .fasta
for i in [".fa",".fas"]:
	if InputGenome.endswith(i):
		os.rename(InputGenome,InputGenome.replace(i,".fasta"))
		InputGenome=InputGenome.replace(i,".fasta")

# function to run revseq (to complement genome, to find antisense hpRNA)
def revseqRun(genome):
	print("Generating complement of "+genome)
	cmd="revseq -sequence "+genome+" -outseq "+genome.replace(".fasta","_comp.fasta")+" -norev"
	subprocess.call(cmd,shell=True)
	print("Complement of sequence "+genome+" written to "+genome.replace(".fasta","_comp.fasta"))
	return(genome.replace(".fasta","_comp.fasta"))

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
		if not line=="\n":
			if line.startswith("         "):
				LeftArm=False
			elif line.startswith("   "):
				temp=line.strip("   ").split(" ")
				if LeftArm==True:
					start=int(temp[0])
					LeftStarts.append(start)
					end=int(temp[2])
					LeftEnds.append(end)
				elif LeftArm==False:
					start=int(temp[0])
					RightStarts.append(start)
					end=int(temp[2])
					RightEnds.append(end)
			else:
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
	BedOutput=""
	for i in range(len(Locations)):
		BedOutput+=Locations[i]+"\t"+str(LeftStarts[i])+"\t"+str(LeftEnds[i])+"\t"+Names[i]+"_left\t"+str(Scores[i])+"\n"+Locations[i]+"\t"+str(RightStarts[i])+"\t"+str(RightEnds[i])+"\t"+Names[i]+"_right\t"+str(Scores[i])+"\n"
	output = open(results.replace(".out",".bed"),"wt")
	output.write(BedOutput)
	output.close()
	print("Bed file of hpRNAs written to "+results.replace(".out",".bed"))
	return()

# master function to: run einverted
def hpRNAfind(input):
#	inputcomp=revseqRun(genome=input)
	einvertedResults=einvertedRun(genome=input)
#	einvertedResultsComp=einvertedRun(genome=inputcomp)
	einvertedParse(results=einvertedResults)
#	einvertedParse(results=einvertedResultsComp)	
	return("hpRNA analysis complete for "+input)

hpRNAfind(input=InputGenome)

