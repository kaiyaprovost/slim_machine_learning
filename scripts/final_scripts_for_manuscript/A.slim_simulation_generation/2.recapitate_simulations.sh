#!/bin/bash

## this script inputs files that look like prefix-TIMESTAMP-i.trees.gz
## and outputs files that look like prefix-TIMESTAMP-i-recap_mu-recomb-Ne-scaling.trees or prefix-TIMESTAMP-i-recap_mu-recomb-Ne-scaling.vcf 
## where mu, recomb, NE, and scaling are given by the values below


newfilename = path+prefix+"-recap_"+str(mu)+"-"+str(recomb)+"-"+str(Ne)+"-"+str(scaling)+"."+suffix
vcffilename = prefix+"-recap_"+str(mu)+"-"+str(recomb)+"-"+str(Ne)+"-"+str(scaling)+".vcf"


export TMPDIR=!/tempfiles/

## default values -- change depending on what you set in previous simulations
mu="1e-7"
recomb="1e-8"
Ne="1000"
scaling="0.02"


for fgz in *trees.gz; do
	gunzip -f $fgz;

	for f in *trees; do

		## these arguments are pulling the file information from the filenames
		## if file names contain or less "-" characters then the information can become shifted
		## if needed, edit the column numbers specified by -f below 
		filename=`echo $f | cut -f1 -d'-'`
		TIMESTAMP=`echo $f | cut -f2 -d'-'`
		suffix=`echo $f | cut -f3 -d'-'`
		i=`echo $suffix | cut -f1 -d'.'`


		echo "${mu} ${recomb} ${Ne}"

		if [ -f "$filename-$TIMESTAMP-$i-recap.trees" ]; then
			echo "$filename-$TIMESTAMP-$i-recap.trees found."
		else
			echo "$filename-$TIMESTAMP-$i-recap.trees not found."
			command="python '''2.1.mutate_and_recapitate.py''' \
			$f \
			$mu \
			$TIMESTAMP \
			$recomb \
			$Ne \
			$scaling ;"
			echo "----------------------------"
			echo "Command is:"
			echo $command
			echo "----------------------------"
			eval $command
			echo "done" 
		fi
	done

	gzip $f
done;

