# This is a Python recipe, to be run after the section 16.1 recipe

import subprocess, msprime, pyslim
import matplotlib.pyplot as plt
import numpy as np
import sys, os, glob

try:
	filename = sys.argv[1]
	print("\tFile is: ",filename)
except:
	print("Filename not given, quitting")
	exit()
	
try:
	mu = float(sys.argv[2])
	print("\tMutation rate (mu, -m) is: ",mu)
except:
	print("Mutation rate not given, quitting")
	exit()

try:
	timestamp = int(sys.argv[3])
	print("\tRandom seed (timestamp) is: ",timestamp)
except:
	print("Random seed/timestamp not given, quitting")
	exit()

splitfile = filename.split("/")

file = splitfile[-1]
path = "/".join(splitfile[:-1])+"/"

if path=="/":
	path=""

prefix = ".".join(file.split(".")[:-1])
suffix = file.split(".")[-1]

newfilename = path+prefix+"-overlaid."+suffix

print("\tLoading")
ts = pyslim.load(filename)
print("\tAdding Mutations")
mutated = msprime.mutate(ts, rate=mu, random_seed=timestamp, keep=True)
print("\tPrinting to file:",newfilename)
mutated.dump(newfilename)

with open(prefix+"-overlaid.vcf", "w") as vcf_file:
	mutated.write_vcf(vcf_file, 2)
