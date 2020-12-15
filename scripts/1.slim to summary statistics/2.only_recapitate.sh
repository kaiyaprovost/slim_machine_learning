#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=99999:00:00
#PBS -N recapitate
#PBS -j oe
#PBS -m ae
#PBS -M kprovost@amnh.org
#PBS -k oe

export TMPDIR=/home/kprovost/nas5/tempfiles/

#source activate py36
cd  "/home/kprovost/nas5/slim_osg/subdirectory/"

for fgz in *trees.gz; do

gunzip -f $fgz;

for f in *trees; do

#for f in migrate-0.1-popsize-200000-mut-2.21e-9-gen-000000-recom-1e-8-ibd-1-seccon-0-pop-2-1582230455-1*trees; do
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
scaling="0.02"

echo "${mu} ${recomb} ${Ne}"

if [ -f "$filename-$TIMESTAMP-$i-recap.trees" ]
then
echo "$filename-$TIMESTAMP-$i-recap.trees found."
else

echo "$filename-$TIMESTAMP-$i-recap.trees not found."

command="python '''/home/kprovost/nas5/slim_osg/Overlaying_neutral_mutations_AND_recapitate_klp_OSG.py''' \
$f \
$mu \
$TIMESTAMP \
$recomb \
$Ne \
$scaling ;"

echo "----------------------------"
echo "Command is:"
echo $command
echo "----------------------------"

eval $command

echo "done" 

fi

done

gzip $f

done;

