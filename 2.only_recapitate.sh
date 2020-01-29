#!/bin/bash

source activate py36
mu="4.42e-8"
recomb="1e-8"
Ne=4000
echo "${mu} ${recomb} ${Ne}"

cd  "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TREES/"

for f in *.trees; do
filename=`echo $f | cut -f1 -d'-'`
TIMESTAMP=`echo $f | cut -f2 -d'-'`
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

mv $filename-$TIMESTAMP-$i-recap.trees* "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TREES/RECAP/"
mv $f* "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TREES/NORECAP/"
mv *png* "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TREES/PNGS/"
mv $filename-$TIMESTAMP-$i-recap.vcf "/home/kprovost/nas2/Analysis_SLiM/FINISHED/VCFS/"

done

echo "

# starting to zip recap"

cd "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TREES/RECAP/"
for i in *trees; do echo "zipping ${i}"; gzip -f $i; done;

echo "

starting to zip non recap"

cd "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TREES/NORECAP/"
for i in *trees; do echo "zipping ${i}"; gzip -f $i; done;
