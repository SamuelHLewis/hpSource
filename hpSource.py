#!/usr/bin/env python3

import argparse
import subprocess
import sys
import os
import re
import shutil
from Bio import SeqIO

# argument parsing
parser=argparse.ArgumentParser(description="Read arguments")
parser.add_argument("-g","--inputgenome",type=str,help="Input genome file")
parser.add_argument("-s","--inputsrna",type=str,help="Input sRNA read file")
parser.add_argument("-c","--cores",type=int,help="Number of cores (default=1)")
parser.add_argument("-p","--gap",type=int,help="Gap penalty (default=12)")
parser.add_argument("-t","--threshold",type=int,help="Minimum score threshold (default=50)")
parser.add_argument("-m","--match",type=int,help="Match score (default=3)")
parser.add_argument("-a","--mismatch",type=int,help="Mismatch score (default=-4)")
parser.add_argument("-f","--fractioncoverage",type=str,help="Threshold for fraction of total coverage contributed by dominant strand (default=0.8)")
parser.add_argument("-d","--distance",type=int,help="Maximum distance (in bases) between hpRNA annotations to be merged into one RNA structure (default=1000")
args=parser.parse_args()
# default values
Cores=1
EinvGap=12
EinvThreshold=50
EinvMatch=3
EinvMismatch=-4
FractionCoverage=0.8
Distance=1000
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
# parse cores
if args.cores is None:
	print("Using default cores ("+str(Cores)+")")
else:
	if args.cores >= 1:
		Cores=int(args.cores)
		print("Cores = "+str(Cores))
	else:
		print("ERROR: cores (-c) must be an integer of 1 or more")
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
# parse fraction of coverage required on dominant strand
if args.fractioncoverage is None:
	print("Using default threshold fraction for dominant strand coverage ("+str(FractionCoverage)+")")
else:
	if 0.5 < float(args.fractioncoverage) <= 1:
		FractionCoverage=float(args.fractioncoverage)
		print("Threshold fraction for dominant strand coverage  = "+str(FractionCoverage))
	else:
		print("ERROR: threshold fraction for dominant strand coverage score (-f) must be above 0.5 and below or equal to 1")
		sys.exit(0)

for i in [".fa",".fas"]:
	if InputGenome.endswith(i):
		os.rename(InputGenome,InputGenome.replace(i,".fasta"))
		InputGenome=InputGenome.replace(i,".fasta")
# parse distance between hpRNA predictions
if args.distance is None:
	print("Using default distance ("+str(Distance)+")")
else:
	if args.distance >= 0:
		Distance=int(args.distance)
		print("Distance = "+str(Distance)+" bases")
	else:
		print("ERROR: distance (-d) must be an integer of 0 or more")
		sys.exit(0)

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
	# NB: initial run is done with very low threshold (minimum 2 adjacent basepairs) because these scores will be modified later to account for G:U wobble pairs, so want to be as permissive as possible at this stage
	cmd="einverted -sequence "+genome+" -gap "+str(EinvGap)+" -threshold "+str(EinvMatch*2)+" -match "+str(EinvMatch)+" -mismatch "+str(EinvMismatch)+" -outfile ./hpSource/"+genome.replace(".fasta",".out")+" -outseq ./hpSource/"+genome.replace(".fasta","_einverted.fasta")
	subprocess.call(cmd,shell=True)
	print("Output of einverted written to ./hpSource/"+genome.replace(".fasta",".out")+" & ./hpSource/"+genome.replace(".fasta","_einverted.fasta"))
	return("./hpSource/"+genome.replace(".fasta",".out"))

# function to calculate score of hairpin, taking into account G:U bases
def einvertedScore(left,right):
	# NB: "g" can pair with "c" or "t", and "t" can pair with "a" or "g", which takes wobble pairs into account
	ScoreDict={
	"a":["t"],
	"c":["g"],
	"g":["c","t"],
	"t":["a","g"],
	"-":["-"]
	}
	if len(left)!=len(right):
		print("ERROR: hairpin arms are different lengths")
		sys.exit(0)
	# calculate score for hpRNA based on the match score and mismatch penalty
	match=0
	mismatch=0
	for base in range(len(left)):
		if right[base] in ScoreDict[left[base]]:
			match+=1
		else:
			mismatch+=1
	score=(match*EinvMatch)+(mismatch*EinvMismatch)
	return(score)

# function to convert einverted output to gff
def einvertedParse(results,reversecomp=False):
	Locations=[]
	Names=[]
	Scores=[]
	LeftStarts=[]
	LeftEnds=[]
	RightStarts=[]
	RightEnds=[]
	# if parsing results from the reverse complement, need to log the total length of each chromosome to convert coordinates back to original (i.e. forward strand) orientation
	if reversecomp:
		CorrectionFactors={}
		for chromosome in SeqIO.parse(InputGenome,"fasta"):
			CorrectionFactors[chromosome.id]=len(chromosome)
			print("Chromosome "+chromosome.id+" is "+str(len(chromosome))+" nt long")
	print("Parsing einverted output file: "+results)
	linecount=0
	hpCount=0
	chrom=""
	leftarm=""
	rightarm=""
	for line in open(results,"r"):
		linecount+=1
		if linecount==1:
			hpCount+=1
		elif linecount==2:
			if chrom==line.split(":")[0]:
				hpCount+=1
			else:
				chrom=line.split(":")[0]
				hpCount=1
		elif linecount==3:
			temp=line.strip("   ").split(" ")
			# if parsing the reverse complement, convert the left start and end coordinates back to the forward orientation (NB: using the "end" coordinate as the start and "start" coordinate as the end for reverse complement)
			if reversecomp:
				leftstart=CorrectionFactors[chrom]-int(temp[2])
				# +1 to end position to account for BED format 0-based indexing
				leftend=CorrectionFactors[chrom]-int(temp[0])+1
			else:
				# -1 from start position to account for BED format 0-based indexing
				leftstart=int(temp[0])-1
				leftend=int(temp[2])
			leftarm=temp[1]
		elif linecount==4:
			None
		elif linecount==5:
			temp=line.strip("   ").split(" ")	
			# if parsing the reverse complement, convert the right start and end coordinates back to the forward orientation (NB: using the "end" coordinate as the start and "start" coordinate as the end for reverse complement)
			if reversecomp:
				rightstart=CorrectionFactors[chrom]-int(temp[0])
				# +1 to end position to account for BED format 0-based indexing
				rightend=CorrectionFactors[chrom]-int(temp[2])+1
			else:
				# -1 from start position to account for BED format 0-based indexing
				rightstart=int(temp[2])-1
				rightend=int(temp[0])
			rightarm=temp[1]
			# calculate the score for this hpRNA using the custom scorer, which takes into account G:U wobble pairs
			score=einvertedScore(left=leftarm,right=rightarm)
			# only log the details of the hpRNA if it passes the threshold
			if score>=EinvThreshold:
				Scores.append(score)
				Locations.append(chrom)
				Names.append(chrom+"_hpRNA"+str(hpCount))
				LeftStarts.append(leftstart)
				LeftEnds.append(leftend)
				RightStarts.append(rightstart)
				RightEnds.append(rightend)
			else:
				None	
			linecount=0
	# output all annotations passing the score threshold in bed format (column order = Chromosome,Start,Stop,Name,Score,Strand)
	BedOutput=""
	if reversecomp:
		for i in range(len(Locations)):
			BedOutput+=Locations[i]+"\t"+str(LeftStarts[i])+"\t"+str(LeftEnds[i])+"\t"+Names[i]+"_left\t"+str(Scores[i])+"\t-\n"+Locations[i]+"\t"+str(RightStarts[i])+"\t"+str(RightEnds[i])+"\t"+Names[i]+"_right\t"+str(Scores[i])+"\t-\n"
	else:
		for i in range(len(Locations)):
			BedOutput+=Locations[i]+"\t"+str(LeftStarts[i])+"\t"+str(LeftEnds[i])+"\t"+Names[i]+"_left\t"+str(Scores[i])+"\t+\n"+Locations[i]+"\t"+str(RightStarts[i])+"\t"+str(RightEnds[i])+"\t"+Names[i]+"_right\t"+str(Scores[i])+"\t+\n"
	output = open(results.replace(".out",".bed"),"wt")
	output.write(BedOutput)
	output.close()
	# append this bed output to a summary file of all hpRNA predictions for both strands
	cmd="cat "+results.replace(".out",".bed")+" >> ./hpSource/premapping_unsorted.bed"
	subprocess.call(cmd,shell=True)
	# sort the bed file and remove the unsorted bedfile
	cmd = "bedtools sort -i ./hpSource/premapping_unsorted.bed > ./hpSource/premapping.bed"
	subprocess.call(cmd,shell=True)
	os.remove("./hpSource/premapping_unsorted.bed")
	print("Bed file of hpRNAs written to "+results.replace(".out",".bed")+" and concatenated to ./hpSource/premapping.bed")
	return()

# function to map sRNAs to a genome
def sRNAmap(srna,genome,cores):
	if srna.endswith(".gz"):
		print("Input file "+srna+" is compressed - uncompressing...")
		cmd="gunzip "+srna
		subprocess.call(cmd,shell=True)
		print("Input file "+srna+" uncompressed")
		srna=srna.replace(".gz","")
	# make bowtie2 database for input genome
	print("Making bowtie2 database for genome "+genome)
	cmd="bowtie2-build "+genome+" ./hpSource/"+genome.replace(".fasta","")
	subprocess.call(cmd,shell=True)
	# map reads to forward strand
	print("Mapping "+srna+" to forward strand of "+genome)
	cmd="bowtie2 --fast --norc -c "+str(cores)+" -q -x ./hpSource/"+genome.replace(".fasta","")+" -U "+srna+" -S ./hpSource/MappedF.sam"
	subprocess.call(cmd,shell=True)
	print("Reads mapped")
	print("Converting SAM to BAM")
	cmd="samtools view -bS -u ./hpSource/MappedF.sam > ./hpSource/MappedF.bam"
	subprocess.call(cmd,shell=True)
	print("SAM converted to BAM")
	os.remove("./hpSource/MappedF.sam")
	# map reads to reverse strand
	print("Mapping "+srna+" to reverse strand of "+genome)
	cmd="bowtie2 --fast --nofw -c "+str(cores)+" -q -x ./hpSource/"+genome.replace(".fasta","")+" -U "+srna+" -S ./hpSource/MappedR.sam"
	subprocess.call(cmd,shell=True)
	print("Reads mapped")
	# convert sam to bam, and remove sam
	print("Converting SAM to BAM")
	cmd="samtools view -bS -u ./hpSource/MappedR.sam > ./hpSource/MappedR.bam"
	subprocess.call(cmd,shell=True)
	print("SAM converted to BAM")
	os.remove("./hpSource/MappedR.sam")
	return()

# function to count reads mapping to each feature in a gff in a strand-specific manner
def CoverageCalculator(bed,bam):
	print("Counting reads in "+bam+" overlapping features in "+bed)
	CoverageFile="./hpSource/"+bed.split("/")[-1].strip(".bed")+"_"+bam.split("/")[-1].strip(".bam")+".counts"
	cmd="bedtools coverage -s -c -a "+bed+" -b "+bam+" > "+CoverageFile
	subprocess.call(cmd,shell=True)
	print("Coverage written to "+CoverageFile)
	return(CoverageFile)

# function to screen out annotations based on the proportion of reads mapping to the dominant strand
def CoverageScreener(forward,reverse):
	print("Screening hpRNAs based on coverage")
	InputLines=[]
	CandidatesF=[]
	CountsF=[]
	CandidatesR=[]
	CountsR=[]
	ValidatedLines=[]
	# read in names of each hpRNA and coverage on either strand
	for line in open(forward,"r"):
		InputLines.append(line.strip("\n"))
		CandidatesF.append(line.split("\t")[3])
		CountsF.append(int(line.split("\t")[-1]))
	for line in open(reverse,"r"):
		CandidatesR.append(line.split("\t")[3])
		CountsR.append(int(line.split("\t")[-1]))
	# work out which strand has higher coverage (dominant strand), calculate what proportion of total coverage this strand accounts for, and retain only those hpRNAs for which this proportion is equal to or more than the specified fraction of coverage contributed by the dominant strand FractionCoverage
	for i in range(len(CandidatesF)):
		if CandidatesF[i]!=CandidatesR[i]:
			print("ERROR: clash of names when parsing forward and reverse counts")
			sys.exit(0)
		else:
			if CountsF[i]>CountsR[i]:
				CoverageProp=CountsF[i]/(CountsF[i]+CountsR[i])
			elif CountsR[i]>CountsF[i]:
				CoverageProp=CountsR[i]/(CountsR[i]+CountsF[i])
			else:
				CoverageProp=0.5
			print("Dominant strand coverage for "+CandidatesF[i]+" = "+str(CoverageProp))
			if CoverageProp>FractionCoverage:
				ValidatedLines.append(InputLines[i])
				print(CandidatesF[i]+" passes validation")
	# output bed file of coverage-validated hpRNA annotations
	ValidatedBed=''
	for i in ValidatedLines:
		temp=i.split("\t")
		ValidatedBed+=temp[0]+"\t"+temp[1]+"\t"+temp[2]+"\t"+temp[3]+"\t"+temp[4]+"\t"+temp[5]+"\n"
	output=open("./hpSource/postmapping_unsorted.bed","wt")
	output.write(ValidatedBed)
	output.close()
	# sort the bed file and remove the unsorted bedfile
	cmd = "bedtools sort -i ./hpSource/postmapping_unsorted.bed > ./hpSource/postmapping.bed"
	subprocess.call(cmd,shell=True)
	os.remove("./hpSource/postmapping_unsorted.bed")
	print("Coverage-screened hpRNA annotations written to ./hpSource/postmapping.bed")
	return()

# function to merge neighbouring hpRNA annotations into one annotation (for RNAfold)
def NeighbourMerger(bedfile,distance=100):
	print("Merging neighbouring hpRNA annotations within "+str(distance)+" bases of each other in file "+bedfile)
	cmd = "bedtools merge -i "+bedfile+" -s -d "+str(distance)+" -c 4,1 -o collapse,count > ./hpSource/merged_wrongorder.bed"
	subprocess.call(cmd,shell=True)
	# the initial ouput of merge is in the wrong order for parsing by genome viewers, so rearrange the columns to put the strand at the end
	cmd="awk \'BEGIN {FS=OFS=\"\t\"} {print $1,$2,$3,$5,$6,$4}\' ./hpSource/merged_wrongorder.bed | sed s/\"\,\"/\"_\"/ > ./hpSource/merged.bed"
	subprocess.call(cmd,shell=True)
	os.remove("./hpSource/merged_wrongorder.bed")
	return()

# function to run RNAfold on a bedfile
def RNAfolder(genome,bedfile):
	os.mkdir("./hpSource/RNAfold/")
	print("Extracting fasta sequences from "+genome+" and "+bedfile)
	cmd="bedtools getfasta -s -fi "+genome+" -bed "+bedfile+" -fo ./hpSource/RNAfold/merged.fasta"
	subprocess.call(cmd,shell=True)
	cmd="RNAfold -i ./hpSource/RNAfold/merged.fasta -o ./hpSource/RNAfold/"
	subprocess.call(cmd,shell=True)
	return()

# master function to find hpRNAs
def hpRNAfind(input):
	if os.path.isdir("hpSource") is True:
		shutil.rmtree("hpSource")
		print("Old hpSource output files removed")
	os.makedirs("hpSource")
	inputRC=revseqRun(genome=input)
	einvertedResults=einvertedRun(genome=input)
	einvertedParse(results=einvertedResults)
	einvertedResultsRC=einvertedRun(genome=inputRC)
	einvertedParse(results=einvertedResultsRC,reversecomp=True)
	sRNAmap(srna=InputsRNA,genome=InputGenome,cores=Cores)
	Fcounts=CoverageCalculator(bed="./hpSource/premapping.bed",bam="./hpSource/MappedF.bam")
	Rcounts=CoverageCalculator(bed="./hpSource/premapping.bed",bam="./hpSource/MappedR.bam")
	CoverageScreener(forward=Fcounts,reverse=Rcounts)
	NeighbourMerger(bedfile="./hpSource/postmapping.bed",distance=Distance)
	RNAfolder(genome=input,bedfile="./hpSource/merged.bed")
	return("hpRNA analysis complete for "+input)

hpRNAfind(input=InputGenome)



