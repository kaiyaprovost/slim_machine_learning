import glob
import sys
import os
import shutil

for file in glob.iglob("/home/kprovost/nas5/slim_osg/subdirectory/*/*bases.vcf_missing*"):
	dirsplit=file.split("/")
	directory="/".join(dirsplit[:-1])
	filename=dirsplit[-1]
	tree=(filename.split("-recap")[0])+".trees.gz"
	treefile=directory+"/"+tree
	newtree=directory+"/DONE/"+tree
	if os.path.exists(treefile):
		shutil.move(treefile,newtree)

for file in glob.iglob("/home/kprovost/nas5/slim_osg/subdirectory/WRINGMIGRATE/*/*bases.vcf_missing*"):
	dirsplit=file.split("/")
	directory="/".join(dirsplit[:-1])
	filename=dirsplit[-1]
	tree=(filename.split("-recap")[0])+".trees.gz"
	treefile=directory+"/"+tree
	newtree=directory+"/DONE/"+tree
	if os.path.exists(treefile):
		shutil.move(treefile,newtree)
