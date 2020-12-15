#!/bin/bash

cd ~/Dropbox\ \(AMNH\)/Dissertation/CHAPTER1_REVIEW/CFB_review_J_Biogeo/snake_fasta_files/Pituophis_vcf/
#for vcffirst in model2*vcf.gz; do 
#vcffile=${vcffirst#./VCFS/}
#vcffile=${vcffirst%.gz}
#gunzip -f $vcffirst
for vcffile in *vcf; do
echo $vcffile
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.fulltemp
numtosample=36 ## atrox 34, getula 26, pituophis 36
#!/bin/bash
if [ -f "$fulltempfile" ]
then
echo "$fulltempfile found."
#mv $vcffile "/home/kprovost/nas2/Analysis_SLiM/FINISHED_SCALED/VCFS/FINISHED/
else
echo "$fulltempfile not found."
perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.pl" $vcffile $fulltempfile $numtosample
fi
gzip -f $vcffile
done	
