#!/bin/bash

source activate py36
cd  "/home/kprovost/nas2/Analysis_SLiM/"

for f in *-?.trees; do
#filename=`echo $f | cut -f1 -d'-'`
#TIMESTAMP=`echo $f | cut -f2 -d'-'`
#suffix=`echo $f | cut -f3 -d'-'`
#i=`echo $suffix | cut -f1 -d'.'`

filename=`echo $f | cut -f1-18 -d'-'`
TIMESTAMP=`echo $f | cut -f19 -d'-'`
suffix=`echo $f | cut -f20 -d'-'`
i=`echo $suffix | cut -f1 -d'.'`

mu=`echo $f | cut -f6-7 -d'-'`
recomb=`echo $f | cut -f11-12 -d'-'`
Ne=`echo $f | cut -f4 -d'-'`
echo "${mu} ${recomb} ${Ne}"

if [ -f "$filename-$TIMESTAMP-$i-recap.trees" ]
then
echo "$filename-$TIMESTAMP-$i-recap.trees found."
else

echo "$filename-$TIMESTAMP-$i-recap.trees not found."

command="python '''/home/kprovost/nas2/Analysis_SLiM/Overlaying_neutral_mutations_AND_recapitate_klp.py''' \
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

gzip -f $filename-$TIMESTAMP-$i-recap*trees
gzip -f $f 

mv $f.gz /home/kprovost/nas2/Analysis_SLiM/TREES/
mv *recap*trees.gz /home/kprovost/nas2/Analysis_SLiM/RECAPTREES/


#mv $filename-$TIMESTAMP-$i-recap.trees* "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TREES/RECAP/"
#mv $f* "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TREES/NORECAP/"
#mv *png* "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TREES/PNGS/"
#mv $filename-$TIMESTAMP-$i-recap.vcf "/home/kprovost/nas2/Analysis_SLiM/FINISHED/VCFS/"

done

#echo "

## starting to zip recap"

#cd "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TREES/RECAP/"
#for i in *trees; do echo "zipping ${i}"; gzip -f $i; done;

#echo "

#starting to zip non recap"

#cd "/home/kprovost/nas2/Analysis_SLiM/FINISHED/TREES/NORECAP/"
#for i in *trees; do echo "zipping ${i}"; gzip -f $i; done;
