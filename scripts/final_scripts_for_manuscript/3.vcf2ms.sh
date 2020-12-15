#!/bin/bash

for vcffile in *vcf; do
echo $vcffile
prefix=${vcffile%.vcf}
fulltempfile=${prefix}.fulltemp

## if you have a different number of samples due to running this on empirical data, change this 
numtosample=20 

#!/bin/bash
if [ -f "$fulltempfile" ]
then
echo "$fulltempfile found."
else
echo "$fulltempfile not found."
perl "vcf2MS.pl" $vcffile $fulltempfile $numtosample
fi
gzip -f $vcffile
done	
