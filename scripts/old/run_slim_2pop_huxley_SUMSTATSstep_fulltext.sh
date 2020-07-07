#!/bin/bash

#TIMESTAMP=`date +"%s"`

#echo "TIMESTAMP:" 
#echo $TIMESTAMP


# /Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/cuvier/slim-test/

cd /Users/kprovost/Documents/Dissertation/CHAPTER1_REVIEW/SLIM/runs/all_done/TEXTFILES/
for d in */*/; do echo $d; bash /Users/kprovost/Documents/Dissertation/CHAPTER1_REVIEW/SLIM/scripts/run_slim_2pop_huxley_SUMSTATSstep_fulltext.sh -d $d; echo "##########"; done
# for d in slim*/TEMPFILES/; do echo $d; bash /Users/kprovost/Documents/Dissertation/CHAPTER1_REVIEW/SLIM/scripts/run_slim_2pop_huxley_SUMSTATSstep.sh -d $d; echo "##########"; done
# for d in slim*/; do echo $d; mv ${d}TEMPFILES/*stats ${d}STATS/; mv ${d}STATS/*sumstats* ${d}; done;

## this script only needs a directory -d (or optionally a directory full of directories) 

while getopts d: OPTION
do
case "$OPTION" in
d)
echo "The value of -d (directory) is $OPTARG"
directory=$OPTARG
#exit
;;	
esac
done

cd $directory

for f in *.txt; do
name=${f%".txt"}
namefull=$name
echo "NAME: ${name}"

## do all of them
echo "(Output to ${namefull}.sumstats.stats)"
cat $f | /Users/kprovost/Documents/Classes/Machine_Learning/popGenMachineLearningExamples-master/msdir/sample_stats >> $namefull.sumstats.stats 
echo " " >> $namefull.sumstats.stats 
echo

done 

