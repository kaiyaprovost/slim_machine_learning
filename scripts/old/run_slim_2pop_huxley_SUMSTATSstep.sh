#!/bin/bash

#TIMESTAMP=`date +"%s"`

#echo "TIMESTAMP:" 
#echo $TIMESTAMP


# /Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/cuvier/slim-test/

# for d in slim*/; do echo $d; bash /Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/run_slim_2pop_huxley_SUMSTATSstep.sh -d $d; echo "##########"; done


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

for f in *.withheader.temp; do
	name=${f%".withheader.temp"}
	namefull=${f%-[0-9]*".withheader.temp"}
	
	## do the one including the iteration
	echo "(Output to ${name}.sumstats.stats)"
	cat $f | /Users/kprovost/Documents/Classes/Machine_Learning/popGenMachineLearningExamples-master/msdir/sample_stats >> ${name}.sumstats.stats 
	echo " " >> $name.sumstats.stats 
	
	## do all of them
	echo "(Output to ${namefull}.sumstats.stats)"
	cat $f | /Users/kprovost/Documents/Classes/Machine_Learning/popGenMachineLearningExamples-master/msdir/sample_stats >> $namefull.sumstats.stats 
	echo " " >> $namefull.sumstats.stats 
	echo
	
done 

