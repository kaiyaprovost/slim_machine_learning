#!/bin/bash

source activate py36

mu="4.42e-8"
recomb="1e-8"
Ne=4000

echo "${mu} ${recomb} ${Ne}"

#cd "/home/kprovost/nas2/Analysis_SLiM/FIN_6K/TREES/MODEL1/"
cd  "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/TREES/"
#cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/all_done/NE1000/TREES/CONVERTED/NOTRECAP/
#cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/

#for f in model*trees; do 
for f in *.trees; do


filename=`echo $f | cut -f1 -d'-'`
TIMESTAMP=`echo $f | cut -f2 -d'-'`
#TIMESTAMP=1
suffix=`echo $f | cut -f3 -d'-'`
i=`echo $suffix | cut -f1 -d'.'`

if [ -f "$filename-$TIMESTAMP-$i-recap.trees" ]
then
echo "$filename-$TIMESTAMP-$i-recap.trees found."
else

echo "$filename-$TIMESTAMP-$i-recap.trees not found."

command="python '''/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/scripts/Overlaying_neutral_mutations_AND_recapitate_klp.py''' \
$f \
$mu \
$TIMESTAMP \
$recomb \
$Ne ;"

## edited to include recapitatation 
# command="python "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/scripts/Overlaying_neutral_mutations_AND_recapitate_klp.py \
# $filename-$TIMESTAMP-$i.trees \
# $mu \
# $TIMESTAMP \
# $recomb \
# $Ne ;"

echo "----------------------------"
echo "Command is:"
echo $command
echo "----------------------------"

eval $command

echo "done" 

fi

echo "
moving"

gzip -f $filename-$TIMESTAMP-$i-recap.trees
gzip -f $f 

mv $filename-$TIMESTAMP-$i-recap.trees* "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/TREES/RECAP/"
mv $f* "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/TREES/NORECAP/"
mv *png* "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/TREES/PNGS/"
mv $filename-$TIMESTAMP-$i-recap.vcf "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/VCFS/"

done

echo "

# starting to zip recap"

cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/TREES/RECAP/"
for i in *trees; do echo "zipping ${i}"; gzip -f $i; done;

echo "

starting to zip non recap"

cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/TREES/NORECAP/"
for i in *trees; do echo "zipping ${i}"; gzip -f $i; done;
