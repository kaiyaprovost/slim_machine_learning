#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=5000:00:00
#PBS -N ms2stats
#PBS -j oe
#PBS -m ae
#PBS -M kprovost@amnh.org
#PBS -k oe

## ON HUXLEY

cd "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/"

for f in SIN*/WINDOWS/*/*ms; do
if [ -f "${f}.sumstats.stats" ]
then
echo "${f}.sumstats.stats found. Skipping."
else
## do the one including the iteration
echo "(Output to ${f}.sumstats.stats)"
cat $f | "/home/kprovost/nas3/sample_stats" > ${f}.sumstats.stats 
echo " " >> $f.sumstats.stats 
fi
done 


## ON MACBOOK

# cd "/Users/kprovost/Dropbox (AMNH)/cardinalis vcf/"
# 
# for f in */*/*ms; do
# 
# if [ -f "${f}.sumstats.stats" ]
# then
# echo "${f}.sumstats.stats found. Skipping."
# else
# ## do the one including the iteration
# echo "(Output to ${f}.sumstats.stats)"
# cat $f | "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/scripts/msdir/sample_stats" > ${f}.sumstats.stats 
# echo " " >> $f.sumstats.stats 
# fi
# 
# done 
