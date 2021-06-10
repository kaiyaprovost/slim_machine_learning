#!/bin/bash

cd "~/simulated_statistics/"

for f in *.withheader.subsettemp; do
	name=${f%".withheader.subsettemp"}
	namefull=${f%-[0-9]*".withheader.subsettemp"}

	if [ -f "${name}.subset.stats" ]; then
		echo "${name}.subset.stats found. Skipping."
	else
		## do the one including the iteration
		echo "(Output to ${name}.subset.stats)"
		cat $f | "~/msdir/sample_stats" > ${name}.subset.stats 
		echo " " >> $name.subset.stats 
	fi

	if [ -f " ${namefull}.subset.sumstats.stats" ]; then
		echo " ${namefull}.subset.sumstats.stats found. Skipping."
	else
		## do all of them
		echo "(Output to ${namefull}.subset.sumstats.stats)"
		cat $f | "!/msdir/sample_stats" >> $namefull.subset.sumstats.stats 
		echo " " >> $namefull.subset.sumstats.stats 
		echo	
	fi	

done 