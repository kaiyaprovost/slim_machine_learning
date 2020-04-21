#!/usr/bin/python3 -u

import sys
import os
import glob
import copy
import random

#vcf="/Users/kprovost/Downloads/test_to_missing_data.vcf"

try:
	vcf = sys.argv[1]
	print("\nRead the following file to add missing:")
	print(vcf)
except:
	print("\nNo file to convert given, quitting")
	sys.exit()
	#toconvert = "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/DXY/textfiles/raw/cri_SON_Dxy_persite_chrfix.txt"

try:
	missing = sys.argv[2]
	print("\nRead the proportion missing to add:")
	print(vcf)
except:
	print("\nNo proportion missing given, defaulting to 0.5")
	missing = 0.5

with open(vcf, "r") as infile:
	lines = infile.readlines()

outstring=vcf+"_missing"+str(missing)+".vcf"

if os.path.exists(outstring) == True:
	print("SKIPPING")
	exit()

## header line is #CHROM
## there are going to be 10 header data before the rest of the data

possiblebases=["A","C","T","G","0","1"]
replacelines=copy.deepcopy(lines)
numinds=None
nummiss=None
ntodraw=None
for i in list(range(len(lines))):
	#print(i)
	thisline = lines[i]
	firstchar = thisline[0]
	if firstchar != "#":
		## do stuff
		#print("A")
		split = thisline.split("\t")
		data = split[9:]
		newdata=copy.deepcopy(data)
		#print("B")
		if(numinds==None):
			numinds=len(data)
		if(nummiss==None):
			nummiss=int(numinds*missing)
		if(ntodraw==None):
			ntodraw=list(range(0,numinds,1))
		#print("C")
		## check if the length of the line matches the expected numinds
		if(len(newdata)!=numinds):
			tempinds=len(newdata)
			thismiss=random.sample(list(range(0,tempinds,1)),int(tempinds*missing))
			print("WARNING: Line "+str(i)+" of "+str(len(lines)-1)+" does not have all data!")
		else:
			thismiss=random.sample(ntodraw,nummiss)
		thismiss.sort()
		for index in thismiss:
			dat2sub = newdata[index]
			for x in possiblebases:
				dat2sub = dat2sub.replace(x,".")
			newdata[index] = dat2sub
		newline = split[:9]+newdata
		newline = "\t".join(newline)
		#print("D")
	else:
		newline = thisline
	replacelines[i] = newline


with open(outstring,"w") as outfile:
	outfile.writelines(replacelines)

