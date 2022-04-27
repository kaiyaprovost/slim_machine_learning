#!/bin/bash
#PBS -l select=1:ncpus=1
#PBS -l walltime=5000:00:00
#PBS -N vcf2ms
#PBS -j oe
#PBS -m ae
#PBS -M kprovost@amnh.org
#PBS -k oe

## Updated:	5 April 2022

## ON HUXLEY

cd "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/VCFS/"

for vcffilegz in SIN*/W*/*/*vcf.gz; do
	prefixgz=${vcffilegz%.vcf.gz}
	fulltempfilegz=${prefixgz}.ms
	#N=25
	#if [ ! -f "$fulltempfilegz" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfilegz}" ]; then
		echo $vcffilegz $prefixgz $fulltempfilegz #$N
		vcffile=${vcffilegz%.gz}
		gunzip -f $vcffilegz
		python "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfilegz #$N
	#else echo "exists"; fi; else echo "exists"; fi
	gzip -f $vcffile
done

for vcffile in SIN*/W*/*/*vcf; do
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.ms
#N=25
#if [ ! -f "$fulltempfile" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfile}" ]; then
python "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfile #$N
#else echo "exists"; fi; else echo "exists"; fi;
gzip -f $vcffile
done


# for vcffilegz in */G*/*/*/*vcf.gz; do
# 	prefixgz=${vcffilegz%.vcf.gz}
# 	fulltempfilegz=${prefixgz}.ms
# 	#N=25
# 	#if [ ! -f "$fulltempfilegz" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfilegz}" ]; then
# 		echo $vcffilegz $prefixgz $fulltempfilegz #$N
# 		vcffile=${vcffilegz%.gz}
# 		gunzip -f $vcffilegz
# 		python "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfilegz #$N
# 	#else echo "exists"; fi; else echo "exists"; fi
# 	gzip -f $vcffile
# done
# 
# for vcffile in */G*/*/*/*vcf; do
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# #N=25
# #if [ ! -f "$fulltempfile" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfile}" ]; then
# python "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfile #$N
# #else echo "exists"; fi; else echo "exists"; fi;
# gzip -f $vcffile
# done



# for vcffilegz in BELLII/W*/G*/*vcf.gz; do
# 	prefixgz=${vcffilegz%.vcf.gz}
# 	fulltempfilegz=${prefixgz}.ms
# 	N=18
# 	if [ ! -f "$fulltempfilegz" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfilegz}" ]; then
# 		echo $vcffilegz $prefixgz $fulltempfilegz $N
# 		vcffile=${vcffilegz%.gz}.vcf
# 		gunzip -f $vcffilegz
# 		perl "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfilegz $N
# 	else echo "exists"; fi; else echo "exists"; fi
# 	gzip -f $vcffile
# done
# for vcffile in BELLII/W*/G*/*vcf; do
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# N=18
# if [ ! -f "$fulltempfile" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfile}" ]; then
# perl "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfile $N
# else echo "exists"; fi; else echo "exists"; fi;
# gzip -f $vcffile
# done
# for vcffilegz in BELLII/G*/*/*/*vcf.gz; do
# 	prefixgz=${vcffilegz%.vcf.gz}
# 	fulltempfilegz=${prefixgz}.ms
# 	N=18
# 	if [ ! -f "$fulltempfilegz" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfilegz}" ]; then
# 		echo $vcffilegz $prefixgz $fulltempfilegz $N
# 		vcffile=${vcffilegz%.gz}.vcf
# 		gunzip -f $vcffilegz
# 		perl "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfilegz $N
# 	else echo "exists"; fi; else echo "exists"; fi
# 	gzip -f $vcffile
# done
# for vcffile in BELLII/G*/*/*/*vcf; do
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# N=18
# if [ ! -f "$fulltempfile" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfile}" ]; then
# perl "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfile $N
# else echo "exists"; fi; else echo "exists"; fi;
# gzip -f $vcffile
# done

# for vcffile in BILINEATA/WINDOWS/*/*vcf; do
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# N=24
# if [ ! -f "$fulltempfile" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfile}" ]; then
# perl "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfile $N
# else echo "exists"; fi; else echo "exists"; fi;
# gzip -f $vcffile
# done
# 
# for vcffile in BRU*/WINDOWS/*/*vcf; do
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# N=22
# if [ ! -f "$fulltempfile" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfile}" ]; then
# perl "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfile $N
# else echo "exists"; fi; else echo "exists"; fi;
# gzip -f $vcffile
# done
# 
# for vcffile in CRI*/WINDOWS/*/*vcf; do
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# N=22
# if [ ! -f "$fulltempfile" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfile}" ]; then
# perl "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfile $N
# else echo "exists"; fi; else echo "exists"; fi;
# gzip -f $vcffile
# done
# 
# for vcffile in CUR*/WINDOWS/*/*vcf; do
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# N=21
# if [ ! -f "$fulltempfile" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfile}" ]; then
# perl "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfile $N
# else echo "exists"; fi; else echo "exists"; fi;
# gzip -f $vcffile
# done
# 
# for vcffile in FLA*/WINDOWS/*/*vcf; do
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# N=22
# if [ ! -f "$fulltempfile" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfile}" ]; then
# perl "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfile $N
# else echo "exists"; fi; else echo "exists"; fi;
# gzip -f $vcffile
# done
# 
# for vcffile in FUS*/WINDOWS/*/*vcf; do
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# N=24
# if [ ! -f "$fulltempfile" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfile}" ]; then
# perl "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfile $N
# else echo "exists"; fi; else echo "exists"; fi;
# gzip -f $vcffile
# done
# 
# for vcffile in MEL*/WINDOWS/*/*vcf; do
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# N=23
# if [ ! -f "$fulltempfile" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfile}" ]; then
# perl "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfile $N
# else echo "exists"; fi; else echo "exists"; fi;
# gzip -f $vcffile
# done
# 
# for vcffile in NIT*/WINDOWS/*/*vcf; do
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# N=20
# if [ ! -f "$fulltempfile" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfile}" ]; then
# perl "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfile $N
# else echo "exists"; fi; else echo "exists"; fi;
# gzip -f $vcffile
# done
# 
# for vcffile in SIN*/WINDOWS/*/*vcf; do
# echo $vcffile
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# N=25
# if [ ! -f "$fulltempfile" ]; then if [ ! -f "/vz-nas1-active/ProcessedGenomicReads/EVERY_PLATE/ANGSD/MSFILES/${fulltempfile}" ]; then
# perl "/home/kprovost/nas3/vcf2MS.py" $vcffile $fulltempfile $N
# else echo "exists"; fi; else echo "exists"; fi;
# gzip -f $vcffile
# done


## ON MACBOOK

# cd "/Users/kprovost/Dropbox (AMNH)/CFB_review_J_Biogeo/"
# 
# for vcffile in *Lamp*vcf; do 
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# ## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
# perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.py" $vcffile $fulltempfile 35
# done
# 
# for vcffile in *caten*vcf; do 
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# ## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
# perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.py" $vcffile $fulltempfile 40
# done
# 
# for vcffile in *scutu*vcf; do 
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# ## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
# perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.py" $vcffile $fulltempfile 36
# done
# 
# for vcffile in *atrox*vcf; do 
# prefix=${vcffile%.vcf}
# fulltempfile=${prefix}.ms
# ## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
# perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.py" $vcffile $fulltempfile 44
# done
# cd "/Users/kprovost/Dropbox (AMNH)/CFB_review_J_Biogeo/scutulatus/g_vcfs/"
# for vcffile in *vcf; do
# prefix=${vcffile%.vcf};
# fulltempfile=${prefix}.ms;
# ## 40 catenefier, 35 lampropeltis, 36 scutulatus, 44 atrox
# if [ ! -f "$fulltempfile" ]; then
# perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.py" $vcffile $fulltempfile 36
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
# perl "/Users/kprovost/Documents/Github/slim_machine_learning/scripts/1.slim to summary statistics/vcf2MS.py" $vcffile $fulltempfile 82
# fi
# done