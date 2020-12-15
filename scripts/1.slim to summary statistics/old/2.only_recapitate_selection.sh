#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=99999:00:00
#PBS -N recapitate
#PBS -j oe
#PBS -m ae
#PBS -M kprovost@amnh.org
#PBS -k oe

source activate py36

cd /home/kprovost/nas5/slim_osg/subdirectory/
cd PANMIXIA
#cd GENEFLOW
#cd ISOLATION
#cd SECONDARYCONTACT

export TMPDIR=/home/kprovost/nas5/tempfiles/

# migrate-0.1-popsize-50000-mut-2.21E-10-gen-21000-recom-1.00E-09-ibd-0-seccon-0-pop-2-growth-1.03-2000-100000-selection-s-2000-0.5-1601587069-50.trees.gz

for fgz in *trees.gz; do
gunzip -f $fgz;

for f in *.trees; do
filename=`echo $f | cut -f1-26 -d'-'`
TIMESTAMP=`echo $f | cut -f27 -d'-'`
suffix=`echo $f | cut -f28 -d'-'`
i=`echo $suffix | cut -f1 -d'.'`
mu=`echo $f | cut -f6-7 -d'-'`
recomb=`echo $f | cut -f11-12 -d'-'`
Ne=`echo $f | cut -f4 -d'-'`
scaling="0.02"

echo "${mu} ${recomb} ${Ne}"
if [ -f "$filename-$TIMESTAMP-$i-recap.trees" ]
then
echo "$filename-$TIMESTAMP-$i-recap.trees found."
gzip -f $f 
mv $f.gz ./DONE/

else
echo "$filename-$TIMESTAMP-$i-recap.trees not found."
command="python3 '''/home/kprovost/nas5/slim_osg/Overlaying_neutral_mutations_AND_recapitate_klp_HUX.py''' \
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

gzip -f $f 
mv $f.gz ./DONE/

fi

done;
done;
