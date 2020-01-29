#!/bin/bash

cd "/home/kprovost/nas2/Analysis_SLiM/FINISHED_SCALED/"
#cd ~/Desktop


#files=(model2*21k*vcf)
#for ((i=${#files[@]}-1; i>=0; i++)); do
#for ((i=${#files[@]}-1; i>=0; i--)); do
#vcffile="${files[$i]}"


for vcffirst in ./VCFS/*vcf; do 

vcffile=${vcffirst#./VCFS/}

echo "#####
starting vcf conversion"

echo $vcffile
## model3_isolation_6k-1551823140-1-overlaid.vcf

prefix=${vcffile%.vcf}
fulltempfile=${prefix}.fulltemp
subsettemp=${prefix}.withheader.subsettemp
subsetvcf=${prefix}.subset.vcf
numtosample=20
locsfile="${prefix}.locs"
echo "LOCS: ${locsfile}"

#!/bin/bash
if [ -f "$fulltempfile" ]
then
echo "$fulltempfile found."
#mv $vcffile "/home/kprovost/nas2/Analysis_SLiM/FINISHED_SCALED/VCFS/FINISHED/
else
echo "$fulltempfile not found."

cd "./VCFS/"
#perl /home/kprovost/nas2/Analysis_SLiM/vcf2MS.pl $vcf $fulltempfile $outSamp
perl "/home/kprovost/nas2/Analysis_SLiM/vcf2MS.pl" $vcffile $fulltempfile $outSamp
cd "/home/kprovost/nas2/Analysis_SLiM/FINISHED_SCALED/"

## now need to subset the fulltempfile so that you get whatever sample size you want on each side of the barrier 
## was doing 10-10 

fi


if [ -f "$subsetvcf" ]
then
echo "$subsetvcf found."
else
echo "$subsetvcf not found."

echo "subsetting"
source activate py36
cd "./VCFS/"

## move locsfile tempoerarily
mv "/home/kprovost/nas2/Analysis_SLiM/FINISHED_SCALED/LOCS/${locsfile}" "/home/kprovost/nas2/Analysis_SLiM/FINISHED_SCALED/VCFS/${locsfile}"
python "/home/kprovost/nas2/Analysis_SLiM/FINISHED_SCALED/select_vcf_individuals.py" $vcffile $numtosample $locsfile
cd "/home/kprovost/nas2/Analysis_SLiM/FINISHED_SCALED/"
mv "/home/kprovost/nas2/Analysis_SLiM/FINISHED_SCALED/VCFS/${locsfile}" "/home/kprovost/nas2/Analysis_SLiM/FINISHED_SCALED/LOCS/${locsfile}"


#if [ -f "$subsettemp" ]
#then
#echo "$subsettemp found."
#else
#echo "$subsettemp not found."

#echo "subsetting"
#source activate py36
#cd "./VCFS/"

## move locsfile tempoerarily
#mv "/home/kprovost/nas2/Analysis_SLiM/FINISHED_SCALED/LOCS/${locsfile}" "/home/kprovost/nas2/Analysis_SLiM/FINISHED_SCALED/VCFS/${locsfile}"
#python "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/scripts/select_ms_individuals.py" $fulltempfile $numtosample $locsfile
#cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/"
#mv "/home/kprovost/nas2/Analysis_SLiM/FINISHED_SCALED/VCFS/${locsfile}" "/home/kprovost/nas2/Analysis_SLiM/FINISHED_SCALED/LOCS/${locsfile}"


echo "
done subsetting"

fi


gzip -f $vcffile
gzip -f $fulltempfile

echo "
moving files"
mv $vcffirst* "./VCFS/DONE/"
mv $fulltempfile "./TEMPS/FULLTEMP/"
mv $subsettemp "./TEMPS/"
mv ./VCFS/*temp ./TEMPS/
mv ./VCFS/*loc* ./LOCS/

done	

echo "

starting to zip vcfs"

cd "./VCFS/DONE/"
for i in *vcf; do echo "zipping ${i}"; gzip -f $i; done;


echo "

starting to zip fulltempfiles"

mv ./VCFS/*temp ./TEMPS/

cd "./VCFS/DONE/"
for i in *temp; do echo "zipping ${i}"; gzip -f $i; done;

cd "./TEMPS/FULLTEMPS/"
for i in *temp; do echo "zipping ${i}"; gzip -f $i; done;