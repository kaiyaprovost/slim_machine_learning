#!/bin/bash

cd "/Users/kprovost/Dropbox (AMNH)/cardinalis vcf/"

for f in */*/*ms; do

if [ -f "${f}.sumstats.stats" ]
then
echo "${f}.sumstats.stats found. Skipping."
else
## do the one including the iteration
echo "(Output to ${f}.sumstats.stats)"
cat $f | "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/scripts/msdir/sample_stats" > ${f}.sumstats.stats 
echo " " >> $f.sumstats.stats 
fi

done 
