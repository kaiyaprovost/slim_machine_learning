#double header fixer python script

import sys
import glob
import os
import gzip

path="~/simulated_statistics/"

os.chdir(path)
listfiles=glob.glob("*.vcf*")

for file in listfiles:
	print(file)
	
	suffix=file[-3:]
	if os.path.exists(path+file+"FIX.vcf"):
		print("skipping")
	elif os.path.exists(path+file+"FIX.vcf.gz"):
		print("skipping")
	else:
		if suffix==".gz":
			with gzip.open(path+file,"rb") as infile:
				lines=infile.readlines()
		else:
			with open(path+file,"r") as infile:
				lines=infile.readlines()	
		if len(lines) <= 12:
			print("PROBLEM -- skipping, too few lines")
			os.rename(path+file,path+file+".SMALL")
		else:
			linelendicti = {}
			prev=None
			breakpoints=[]
			lenbreaks=[]
			differences=[0]
			for i in range(len(lines)):
				line = lines[i]
				if suffix==".gz":
					thislength=len(line.split(b"\t"))
				else:
					thislength=len(line.split("\t"))
				entry = linelendicti.get(thislength,None)
				if entry==None:
					linelendicti[thislength] = 1
				else:
					linelendicti[thislength] += 1
				if prev != thislength or i == len(lines)-1:
					if i != 0:
						differences.append(i-breakpoints[-1])
					breakpoints.append(i)
					lenbreaks.append(thislength)
					prev = thislength
			end=differences.index(max(differences))
			headerstarts = [i for i, value in enumerate(lenbreaks) if value == 1] 
			endline = breakpoints[end]+1
			for j in reversed(headerstarts):
				if j < end:
					start=j
					break
			startline=breakpoints[start]
			newlines=lines[startline:endline]
			if suffix==".gz":
				with gzip.open(path+file+"FIX.vcf.gz","wb") as outfile:
					outfile.writelines(newlines)
			else:
				with open(path+file+"FIX.vcf","w") as outfile:
					outfile.writelines(newlines)

print("done")













