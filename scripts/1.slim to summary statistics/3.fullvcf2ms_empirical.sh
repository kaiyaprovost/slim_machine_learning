#!/bin/bash

cd "/Users/kprovost/Dropbox (AMNH)/CFB_review_J_Biogeo/"

for vcffile in *Lamp*vcf; do 
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.pl" $vcffile $fulltempfile 35
done

for vcffile in *caten*vcf; do 
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.pl" $vcffile $fulltempfile 40
done

for vcffile in *scutu*vcf; do 
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.pl" $vcffile $fulltempfile 36
done

for vcffile in *atrox*vcf; do 
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.pl" $vcffile $fulltempfile 44
done
cd "/Users/kprovost/Dropbox (AMNH)/CFB_review_J_Biogeo/scutulatus/g_vcfs/"
for vcffile in *vcf; do
prefix=${vcffile%.vcf};
fulltempfile=${prefix}.ms;
## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
if [ ! -f "$fulltempfile" ]; then
perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.pl" $vcffile $fulltempfile 36
fi
done

cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/CFB_review_J_Biogeo/cardinalis vcf/"
for vcffile in *temp*thr*.vcf; do
echo $vcffile
prefix=${vcffile%.vcf};
fulltempfile=${prefix}.ms;
## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
if [ ! -f "$fulltempfile" ]; then
perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.pl" $vcffile $fulltempfile 82
fi
done