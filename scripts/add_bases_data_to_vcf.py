#!/usr/bin/python3 -u

import sys
import os
import glob
import copy
import random

#vcf="/Users/kprovost/Downloads/TEST.vcf"

try:
	vcf = sys.argv[1]
	print("\nRead the following file to add missing:")
	print(vcf)
except:
	print("\nNo file to convert given, quitting")
	sys.exit()
	#toconvert = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/raw/cri_SON_Dxy_persite_chrfix.txt"

with open(vcf, "r") as infile:
	lines = infile.readlines()

outstring=vcf+"_bases.vcf"

## header line is #CHROM
## there are going to be 10 header data before the rest of the data
## the columns we need are REF and ALT which are cols 4 and 5 (1 index)

possiblebases=["A","C","T","G"]
replacelines=copy.deepcopy(lines)

for i in list(range(len(lines))):
	#print(i)
	thisline = lines[i]
	firstchar = thisline[0]
	if firstchar != "#":
		## do stuff
		split = thisline.split("\t")
		newsplit = copy.deepcopy(split)
		newsplit[3:5]=random.sample(possiblebases,2)
		newline = "\t".join(newsplit)
	else:
		newline = thisline
	replacelines[i] = newline

with open(outstring,"w") as outfile:
	outfile.writelines(replacelines)


