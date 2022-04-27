#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=5000:00:00
#PBS -N ms2stats
#PBS -j oe
#PBS -m ae
#PBS -M kprovost@amnh.org
#PBS -k oe

## Updated:	4 April 2022


## ON HUXLEY

cd "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/"

for fgz in /vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/*/*/*/*ms.gz; do
f=${fgz%.gz}
if [ -f "${f}.sumstats.stats" ]
then
echo "${f}.sumstats.stats found. Skipping."
else
## do the one including the iteration
echo "(Output to ${f}.sumstats.stats)"
gunzip -f $fgz
cat $f | "/home/kprovost/nas3/msdir/sample_stats" > ${f}.sumstats.stats 
echo " " >> $f.sumstats.stats 
gzip -f $f
fi
done 

for f in /vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/*/*/*/*ms; do
if [ -f "${f}.sumstats.stats" ]
then
echo "${f}.sumstats.stats found. Skipping."
else
## do the one including the iteration
echo "(Output to ${f}.sumstats.stats)"
cat $f | "/home/kprovost/nas3/msdir/sample_stats" > ${f}.sumstats.stats 
echo " " >> $f.sumstats.stats 
gzip -f $f
fi
done 

for fgz in /vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/*/*/*ms.gz; do
f=${fgz%.gz}
if [ -f "${f}.sumstats.stats" ]
then
echo "${f}.sumstats.stats found. Skipping."
else
## do the one including the iteration
echo "(Output to ${f}.sumstats.stats)"
gunzip -f $fgz
cat $f | "/home/kprovost/nas3/msdir/sample_stats" > ${f}.sumstats.stats 
echo " " >> $f.sumstats.stats 
gzip -f $f
fi
done

for f in /vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/*/*/*ms; do
if [ -f "${f}.sumstats.stats" ]
then
echo "${f}.sumstats.stats found. Skipping."
else
## do the one including the iteration
echo "(Output to ${f}.sumstats.stats)"
cat $f | "/home/kprovost/nas3/msdir/sample_stats" > ${f}.sumstats.stats 
echo " " >> $f.sumstats.stats 
gzip -f $f
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
# cat $f | "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/sBELpts/msdir/sample_stats" > ${f}.sumstats.stats 
# echo " " >> $f.sumstats.stats 
# fi
# 
# done 
