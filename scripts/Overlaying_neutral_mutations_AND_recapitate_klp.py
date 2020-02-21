## this is a python script adapted from
## recipe 16.2 and 16.10 from the
## SLIM manual

import subprocess, msprime, pyslim
import matplotlib.pyplot as plt
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
	mu=mu*scaling
	recomb=recomb*scaling
	Ne=int(Ne*scaling)

print("\tLoading: "+filename)
ts = pyslim.load(filename) ## do not simplify 

## add the tree_heights function from slim manual 
print("\tDefining Tree Heights")
def tree_heights(ts):
    heights = np.zeros(ts.num_trees + 1)
    for tree in ts.trees():
        if tree.num_roots > 1:  # not fully coalesced
            heights[tree.index] = ts.slim_generation
        else:
            children = tree.children(tree.root)
            real_root = tree.root if len(children) > 1 else children[0]
            heights[tree.index] = tree.time(real_root)
    heights[-1] = heights[-2]  # repeat the last entry for plotting with step
    return(heights)

## plot tree heights before recapacitation
#print("\tPlotting tree heights (before)")
#breakpoints = list(ts.breakpoints())
#heights = tree_heights(ts)
#plt.step(breakpoints, heights, where='post')
#plt.savefig("Before_Recap_"+filename+".png")

# Recapitate!
print("\tRecapitating")
recap = ts.recapitate(recombination_rate=recomb, Ne=Ne, random_seed=timestamp)
## do not output the recapacitated tree 
#recap.dump("recipe_16.10_recap.trees")

# Plot the tree heights after recapitation
# print("\tPlotting tree heights (after): "+"After_Recap_"+filename+".png")
# breakpoints = list(recap.breakpoints())
# heights = tree_heights(recap)
# plt.step(breakpoints, heights, where='post')
# plt.savefig("After_Recap_"+filename+".png")

print("\tAdding Mutations")
mutated = msprime.mutate(recap, rate=mu, random_seed=timestamp, keep=True, model=msprime.InfiniteSites(alphabet=1)) ## nucleotides
## wont do nucleotides -- kelleher and ralph to fix -- ralph receptive to questions with python side 
## can't round positions 

print("\tPrinting to file:",newfilename)
mutated.dump(newfilename)

print("\tCreating VCF")
with open(vcffilename, "w") as vcf_file:
    mutated.write_vcf(vcf_file, 2)
    
    # Plot the tree heights after recapitation
#     print("\tPlotting tree heights (mutated)")
#     breakpoints = list(recap.breakpoints())
#     heights = tree_heights(recap)
#     plt.step(breakpoints, heights, where='post')
#     plt.savefig("After_Mutate_"+filename+".png")
