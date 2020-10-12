#!/bin/bash

#TIMESTAMP=`date +"%s"`

#echo "TIMESTAMP:" 
#echo $TIMESTAMP


cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/TEMPS/"

# while getopts d: OPTION
# do
# 	case "$OPTION" in
# 		d)
# 			echo "The value of -d (directory) is $OPTARG"
# 			directory=$OPTARG
# 			#exit
# 			;;	
# 	esac
# done
# 
# cd $directory

for f in *.withheader.subsettemp; do
name=${f%".withheader.subsettemp"}
namefull=${f%-[0-9]*".withheader.subsettemp"}

if [ -f "${name}.subset.stats" ]
then
echo "${name}.subset.stats found. Skipping."
#mv $f "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/TEMPS/SUBSET/
else
## do the one including the iteration
echo "(Output to ${name}.subset.stats)"
cat $f | "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/scripts/msdir/sample_stats" > ${name}.subset.stats 
echo " " >> $name.subset.stats 
fi

if [ -f " ${namefull}.subset.sumstats.stats" ]
then
echo " ${namefull}.subset.sumstats.stats found. Skipping."
#mv $f "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/TEMPS/SUBSET/
else
## do all of them
echo "(Output to ${namefull}.subset.sumstats.stats)"
cat $f | "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/scripts/msdir/sample_stats" >> $namefull.subset.sumstats.stats 
echo " " >> $namefull.subset.sumstats.stats 
echo	
fi	

#gzip -f ${namefull}.subset.sumstats.stats

#mv $f "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/TEMPS/SUBSET/"
#mv ${name}.subset.stats "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS/SUMSTAT/"
#mv ${namefull}.subset.sumstats.stats* "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS/FULLSUMSTAT/"

done 

#echo "
#
#starting to zip fullsumstat"
#
#cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/STATS/FULLSUMSTAT/"
#for i in *stats; do echo "zipping ${i}"; gzip -f $i; done;