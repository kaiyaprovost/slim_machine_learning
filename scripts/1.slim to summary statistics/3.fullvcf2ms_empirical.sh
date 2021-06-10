#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=5000:00:00
#PBS -N vcf2ms
#PBS -j oe
#PBS -m ae
#PBS -M kprovost@amnh.org
#PBS -k oe

## ON HUXLEY

cd "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/"

for vcffile in BELLII/WINDOWS/*/*vcf; do
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
N=18
if [ ! -f "$fulltempfile" ]; then
perl "/home/kprovost/nas3/vcf2MS.pl" $vcffile $fulltempfile $N
fi
gzip -f $vcffile
done

for vcffile in BILINEATA/WINDOWS/*/*vcf; do
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
N=24
if [ ! -f "$fulltempfile" ]; then
perl "/home/kprovost/nas3/vcf2MS.pl" $vcffile $fulltempfile $N
fi
gzip -f $vcffile
done

for vcffile in BRU*/WINDOWS/*/*vcf; do
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
N=22
if [ ! -f "$fulltempfile" ]; then
perl "/home/kprovost/nas3/vcf2MS.pl" $vcffile $fulltempfile $N
fi
gzip -f $vcffile
done

for vcffile in CRI*/WINDOWS/*/*vcf; do
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
N=22
if [ ! -f "$fulltempfile" ]; then
perl "/home/kprovost/nas3/vcf2MS.pl" $vcffile $fulltempfile $N
fi
gzip -f $vcffile
done

for vcffile in CUR*/WINDOWS/*/*vcf; do
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
N=21
if [ ! -f "$fulltempfile" ]; then
perl "/home/kprovost/nas3/vcf2MS.pl" $vcffile $fulltempfile $N
fi
gzip -f $vcffile
done

for vcffile in FLA*/WINDOWS/*/*vcf; do
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
N=22
if [ ! -f "$fulltempfile" ]; then
perl "/home/kprovost/nas3/vcf2MS.pl" $vcffile $fulltempfile $N
fi
gzip -f $vcffile
done

for vcffile in FUS*/WINDOWS/*/*vcf; do
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
N=24
if [ ! -f "$fulltempfile" ]; then
perl "/home/kprovost/nas3/vcf2MS.pl" $vcffile $fulltempfile $N
fi
gzip -f $vcffile
done

for vcffile in MEL*/WINDOWS/*/*vcf; do
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
N=23
if [ ! -f "$fulltempfile" ]; then
perl "/home/kprovost/nas3/vcf2MS.pl" $vcffile $fulltempfile $N
fi
gzip -f $vcffile
done

for vcffile in NIT*/WINDOWS/*/*vcf; do
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
N=20
if [ ! -f "$fulltempfile" ]; then
perl "/home/kprovost/nas3/vcf2MS.pl" $vcffile $fulltempfile $N
fi
gzip -f $vcffile
done

for vcffile in SIN*/WINDOWS/*/*vcf; do
echo $vcffile
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
N=25
if [ ! -f "$fulltempfile" ]; then
perl "/home/kprovost/nas3/vcf2MS.pl" $vcffile $fulltempfile $N
fi
gzip -f $vcffile
done


## ON MACBOOK

# cd "/Users/kprovost/Dropbox (AMNH)/CFB_review_J_Biogeo/"
# 
# for vcffile in *Lamp*vcf; do 
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# ## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
# perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.pl" $vcffile $fulltempfile 35
# done
# 
# for vcffile in *caten*vcf; do 
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# ## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
# perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.pl" $vcffile $fulltempfile 40
# done
# 
# for vcffile in *scutu*vcf; do 
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# ## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
# perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.pl" $vcffile $fulltempfile 36
# done
# 
# for vcffile in *atrox*vcf; do 
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# ## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
# perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.pl" $vcffile $fulltempfile 44
# done
# cd "/Users/kprovost/Dropbox (AMNH)/CFB_review_J_Biogeo/scutulatus/g_vcfs/"
# for vcffile in *vcf; do
# prefix=${vcffile%.vcf};
# fulltempfile=${prefix}.ms;
# ## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
# if [ ! -f "$fulltempfile" ]; then
# perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.pl" $vcffile $fulltempfile 36
# fi
# done
# 
# cd "/Users/kprovost/Dropbox (AMNH)/Dissertation/CHAPTER1_REVIEW/CFB_review_J_Biogeo/cardinalis vcf/"
# for vcffile in *temp*thr*.vcf; do
# echo $vcffile
# prefix=${vcffile%.vcf};
# fulltempfile=${prefix}.ms;
# ## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
# if [ ! -f "$fulltempfile" ]; then
# perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.pl" $vcffile $fulltempfile 82
# fi
# done