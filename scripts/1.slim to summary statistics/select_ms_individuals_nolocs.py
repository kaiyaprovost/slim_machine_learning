# python 

## usage: python3 select_ms_individuals_nolocs.py $fulltempfile $numtosample

# cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/fulltemps/"
# cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/all_done/actual_scaled_runs"
# for fulltempfile in *fulltemp; do 
# python3 "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/select_ms_individuals_nolocs.py" $fulltempfile 20
# done
# for fulltempfilegz in *fulltemp.gz; do
# gunzip -f $fulltempfilegz
# fulltempfile=${fulltempfilegz%.gz}
# python3 "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/select_ms_individuals_nolocs.py" $fulltempfile 20
# gzip -f $fulltempfile
# done

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
	fulltempfile = "/Users/kprovost/Documents/Dissertation/CHAPTER1_REVIEW/SLIM/runs/all_done/VCFS/6000/ISO/atest/model3_isolation_6k-1551823140-1.fulltemp"

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