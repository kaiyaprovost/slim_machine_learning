# python 

## this script will subset a large number of individuals in an MS file to a small number of individuals 
## usage: python3 select_ms_individuals_nolocs.py $fulltempfile $numtosample

import random 
import numpy as np
import sys

## input the file

try:
	fulltempfile = sys.argv[1]
	print("\nRead the following filename:")
	print(fulltempfile)
except:
	print("\nNo filename given, quitting")
	sys.exit()

try:
	numtosample = int(sys.argv[2])
	print("\nRead the following number to sample:")
	print(numtosample)
except:
	print("\nNo/improper sample number given, defaulting to 20")
	numtosample=20
	
## read in the file 
with open(fulltempfile,"r") as filename:
	lines = filename.readlines()

print("Numlines of fulltempfile: ",len(lines))

## get the data and the metadata
## three lines of metadata
separator = lines[0]
segsites = lines[1]
positions = lines[2]
## actual data
data = lines[3:]

## input the number to sample total and then per pop
## this is num individuals, num genomes is double this
numtosampleperpop = numtosample//2

## get the number of genomes and individuals
numgenomes = len(data)
numinds = numgenomes//2

## set the number of populations -- here always two
numpops = 2 

## extract individuals per populations
numperpop = numinds//2 

## extract individuals per population
firstpop = data[:numperpop*2]
secondpop = data[numperpop*2:]

## randomly select the sample of individuals
topickpop1 = np.random.choice(a=numperpop,size=numtosampleperpop,replace=False)
topickpop2 = numperpop+np.random.choice(a=numperpop,size=numtosampleperpop,replace=False)

"SUBSETTING GENOMES"
## convert individuals to genomes 
genomespop1_fwd = topickpop1*2
genomespop1_rev = 1+topickpop1*2

genomespop2_fwd = topickpop2*2
genomespop2_rev = 1+topickpop2*2

## glue the genome lists together 
genomespop1 = (np.ndarray.tolist(genomespop1_fwd)+np.ndarray.tolist(genomespop1_rev))
genomespop1.sort()

genomespop2 = (np.ndarray.tolist(genomespop2_fwd)+np.ndarray.tolist(genomespop2_rev))
genomespop2.sort()

genomes_all = genomespop1+genomespop2
genomes_all.sort()

## subset the data by those genomes
subsetdata = [data[i] for i in genomes_all]

header = "slim 100 "+str(numtosample)+"\n"

filenameprefix = fulltempfile[:-9]
subsetfilename = filenameprefix+".withheader.subsettemp"

print("WRITING GENOMES")
## write genomes then write locs
with open(subsetfilename,"w") as outfile:
	outfile.write(header)
	outfile.write(separator)
	outfile.write(segsites)
	outfile.write(positions)
	for i in subsetdata: 
		outfile.write(i)