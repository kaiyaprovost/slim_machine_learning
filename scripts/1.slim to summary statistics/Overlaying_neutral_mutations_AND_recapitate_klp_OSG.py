## this is a python script adapted from
## recipe 16.2 and 16.10 from the
## SLIM manual

import subprocess, msprime, pyslim
import numpy as np
import sys, os, glob

try:
    filename = sys.argv[1]
    print("\tFile is: ",filename)
except:
    print("Filename not given, quitting")
    ## model1_panmixia_6k-1557862508-15.trees
    exit()

try:
    mu = float(sys.argv[2])
    print("\tMutation rate (mu, -m) is: ",mu)
except:
    print("Mutation rate not given, quitting")
    ## "1e-7"
    exit()

try:
    timestamp = int(sys.argv[3])
    print("\tRandom seed (timestamp) is: ",timestamp)
except:
    print("Random seed/timestamp not given, quitting")
    ## 1557862508
    exit()

try:
    recomb = float(sys.argv[4])
    print("\tRecombination rate (recomb, -r) is: ",recomb)
except:
    print("Recombination rate not given, quitting")
    ## "1e-8"
    exit()

try:
    Ne = int(sys.argv[5])
    print("\tEff. pop size (Ne, -N) is: ",Ne)
except:
    print("Ne not given, quitting")
    ## 1000 
    exit()

try:
    scaling = float(sys.argv[6])
    print("\tScaling factor (-X) is: ",scaling)
except:
    print("Scaling factor not given, defaulting to 1")
    scaling=1.0

splitfile = filename.split("/")

file = splitfile[-1]
path = "/".join(splitfile[:-1])+"/"

if path=="/":
    path=""

prefix = ".".join(file.split(".")[:-1])
suffix = file.split(".")[-1]

newfilename = path+prefix+"-recap_"+str(mu)+"-"+str(recomb)+"-"+str(Ne)+"-"+str(scaling)+"."+suffix
vcffilename = prefix+"-recap_"+str(mu)+"-"+str(recomb)+"-"+str(Ne)+"-"+str(scaling)+".vcf"

## import os
exists = os.path.isfile(newfilename)

if exists:
    print("\t~~~~~Recap trees exists -- exiting")
    exit()

if(scaling != 1):
    print("\tScaling parameters")
    mu=mu/scaling
    recomb=recomb/scaling
    Ne=int(Ne*scaling)

print("\tLoading: "+filename)
ts = pyslim.load(filename) ## do not simplify 

## check Ne
if Ne!=ts.num_individuals:
    print("ADJUSTING NE TO:",str(ts.num_individuals))
    Ne = int(ts.num_individuals)

# Recapitate!
print("\tRecapitating")
recap = ts.recapitate(recombination_rate=recomb, Ne=Ne, random_seed=timestamp)

print("\tAdding Mutations")
mutated = msprime.mutate(recap, rate=mu, random_seed=timestamp, keep=True, model=msprime.InfiniteSites(alphabet=1)) ## nucleotides

print("\tPrinting to file:",newfilename)
mutated.dump(newfilename)

print("\tCreating VCF")
with open(vcffilename, "w") as vcf_file:
    try:
        #mutated.write_vcf(vcf_file, ploidy=2)
        mutated.write_vcf(vcf_file)    
        #mutated.write_vcf(mutated.individuals_alive_at(0))
    except:
        print("\tSimplifying")
        simple=mutated.simplify()
        simple.write_vcf(vcf_file)



