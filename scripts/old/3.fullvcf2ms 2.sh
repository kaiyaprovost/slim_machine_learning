#!/bin/bash

cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER2_GENOMES/ANALYSIS/called_geno/"
#cd ~/Desktop


#files=(model2*21k*vcf)
#for ((i=${#files[@]}-1; i>=0; i++)); do
#for ((i=${#files[@]}-1; i>=0; i--)); do
#vcffile="${files[$i]}"

for vcffirst in */*vcf; do 

vcffile=${vcffirst#./runs/VCFS/}

echo "#####
starting vcf conversion"

echo $vcffile
## model3_isolation_6k-1551823140-1-overlaid.vcf

prefix=${vcffile%.vcf}
fulltempfile=${prefix}.fulltemp
# subsettemp=${prefix}.withheader.subsettemp
outSamp=18
numtosample=18
# locsfile="${prefix}.locs"
# echo "LOCS: ${locsfile}"

#!/bin/bash
if [ -f "$fulltempfile" ]
then
echo "$fulltempfile found."
#mv $vcffile "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/VCFS/FINISHED/
else
echo "$fulltempfile not found."

#cd "./runs/VCFS/"
#perl /home/kprovost/nas2/Analysis_SLiM/vcf2MS.pl $vcf $fulltempfile $outSamp
perl "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/scripts/vcf2MS.pl" $vcffile $fulltempfile $outSamp
#cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/"
#cd "/Users/kprovost/Documents/folder_for_popgenome/"


## now need to subset the fulltempfile so that you get whatever sample size you want on each side of the barrier 
## was doing 10-10 

fi

# if [ -f "$subsettemp" ]
# then
# echo "$subsettemp found."
# else
# echo "$subsettemp not found."
# 
# echo "subsetting"
# source activate py36
# cd "./runs/VCFS/"
# 
# ## move locsfile tempoerarily
# mv "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/LOCS/${locsfile}" "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/VCFS/${locsfile}"
# python "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/scripts/select_ms_individuals.py" $fulltempfile $numtosample $locsfile
# cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/"
# mv "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/VCFS/${locsfile}" "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/SLIM/runs/LOCS/${locsfile}"
# 
# 
# echo "
# done subsetting"
# 
# fi/


# gzip -f $vcffile
# gzip -f $fulltempfile

# echo "
# moving files"
# mv $vcffirst* "./runs/VCFS/DONE/"
# mv $fulltempfile "./runs/TEMPS/FULLTEMP/"
# mv $subsettemp "./runs/TEMPS/"
# mv ./runs/VCFS/*temp ./runs/TEMPS/
# mv ./runs/VCFS/*loc* ./runs/LOCS/

done	

# echo "
# 
# starting to zip vcfs"
# 
# cd "./runs/VCFS/DONE/"
# for i in *vcf; do echo "zipping ${i}"; gzip -f $i; done;
# 
# 
# echo "
# 
# starting to zip fulltempfiles"
# 
# mv ./runs/VCFS/*temp ./runs/TEMPS/
# 
# cd "./runs/VCFS/DONE/"
# for i in *temp; do echo "zipping ${i}"; gzip -f $i; done;
# 
# cd "./runs/TEMPS/FULLTEMPS/"
# for i in *temp; do echo "zipping ${i}"; gzip -f $i; done;