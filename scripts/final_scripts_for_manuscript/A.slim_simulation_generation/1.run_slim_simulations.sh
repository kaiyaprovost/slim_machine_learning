#!/bin/bash

## ONE OF THESE FOR EACH SIMULATION TYPE 

## example run -- this example will run for 6000 generations with panmixia and no ibd
# bash 1.run_slim_simulations.sh \
# -a 0 \
# -B 1 \
# -c 1 \
# -e "vireotest_slim_TRANSITION_57x87_coarse.txt" \
# -f 1000 \
# -G 1000 \
# -g 1.00 \
# -H 1000 \
# -h "n" \
# -I "0.run_slim_simulations.slim" \
# -i 0 \
# -l 1 \
# -M 0.0 \
# -m "1e-7" \
# -N 1000 \
# -o "panmixia_noibd_6000" \
# -O 6000 \
# -P "" \
# -p slim \
# -r "1e-8" \
# -s 20 \
# -S 100 \
# -t 1 \
# -v 0 \
# -V 5000 \
# -w "0.0" \
# -x 57 \
# -X 0.02 \
# -y 87 
## this example file will produce files with filename resembling panmixia_noibd_6000-TIMESTAMP-i.trees
## where TIMESTAMP is the system time in seconds that the run is initialized 
## and i is the specific run number, ranging from 1 to the value given by -S
## it can also produce a .locs file, .txt file, .log file, and .temp file 

export TMPDIR=~/tempfiles/

TIMESTAMP=`date +"%s"`

echo "TIMESTAMP:" 
echo $TIMESTAMP

## defining arguments					## suggested values
# -a flag # use all pops or not 		# 0 or 1
# -B flag # number of subpops 			# 1 or 0
# -c flag # competition yes/no			# 1 or 0; only works if spatial 1
# -e flag # ENM filename 				# vireotest_slim_27x45.asc
# -f flag # ending pop size 			# 1000
# -G flag # growth generation start		# 1000
# -g flag # growth/contract rate 		# 1 (no population size change), less than 1 (contraction), greater than 1 (growth)
# -H flag # sweep generation start		# 1000
# -h flag # what type of sweep wanted 	# "n" for no sweep, "h" for hard sweep, "s" for soft sweep
# -i flag # ibd yes/no					# 1 or 0; only works if spatial 1
# -I flag # which slim file 			# "0.run_slim_simulations.slim"
# -l flag # spatial yes/no 				# 0 or 1
# -M flag # migration rate 				# 0.1; will accept ranges 0.0 to 1.0, but do not recommend higher than 0.1
# -m flag # mutation rate 				# 1e-7
# -N flag # beginning pop size 			# 1000
# -o flag # name of output file 		# "modeltest"
# -O flag # output generation 			# 2000
# -P flag # path 						# ""; set to path where analyses should be output  
# -p flag # where slim is located 		# "slim"; set to path where slim is found if not found in same folder as this script
# -r flag # recombination rate 			# 1e-8
# -s flag # number genomes to output 	# 20
# -S flag # number of reps to run 		# 100; use 2 for testing
# -t flag # use treeseq or not 			# 1 or 0
# -v flag # secondary contact or not 	# 0 or 1
# -V flag # secondary contact start		# 1000
# -w flag # selection coefficient		# 0 (neutral), less than 0 (deleterious), greater than 0 (beneficial)
# -x flag # niche model X dimension		# 27; do not change if using default ascii file
# -X flag # scaling factor				# 0.02
# -y flag # niche model Y dimension		# 45; do not change if using default ascii file

while getopts p:S:I:P:o:m:r:N:g:f:s:w:h:x:y:e:c:i:l:t:B:a:G:M:O:v:V:H:X: OPTION
do
	case "$OPTION" in
		p)
			echo "The value of -p (path to Slim) is $OPTARG";
			slimPath=$OPTARG;
			#exit
			;;	
		S)
			echo "The value of -S (number of sims to run) is $OPTARG";
			simulations=$OPTARG;
			#exit
			;;
		I)
			echo "The value of -I (SLiM file called) is $OPTARG";
			script=$OPTARG;
			#exit
			;;
		P)
			echo "The value of -P (slim file path) is $OPTARG";
			path=$OPTARG;
			#exit
			;;
		o)
			echo "The value of -o (output filename(s)) is $OPTARG";
			filename=$OPTARG;
			#exit
			;;
		m)
			echo "The value of -m (mu) is $OPTARG";
			mu=$OPTARG;
			#exit
			;;			
		r)
			echo "The value of -r (recombination rate) is $OPTARG";
			recomb=$OPTARG;
			#exit
			;;
		N)
			echo "The value of -N (N individuals in population 1) is $OPTARG";
			N=$OPTARG;
			N_dif=$N #final pop to reach after growth or contract # consider include check that if g_rate under 1, N_dif must be less than N # if g_rate = 1, N_dif not needed?
			#exit
			;;
		g)
			echo "The value of -g (growth or contraction rate) is $OPTARG";
			g_rate=$OPTARG;
			#exit
			;;		
		f)
			echo "The value of -f (final population after growth/contraction) is $OPTARG";
			N_dif=$OPTARG;
			#exit
			;;
		s)
			echo "The value of -s (number of individuals outputted from SLiM) is $OPTARG";
			echo "~~~Note: each individual has 2 genomes output for any given number -- slim multiplies this by 2"
			echo "~~~However: if treeSeq is on, it is not multiplied by 2"
			outSamp=$OPTARG;
			#exit
			;;		
		w)
			echo "The value of -w (selection strength of mutation) is $OPTARG";
			selStrength=$OPTARG;
			#exit
			;;	
		h)
			echo "The value of -h (sweep type) is $OPTARG";
			sweepType=$OPTARG;
			;;
		x)
			echo "The value of -x (spatial x) is $OPTARG";
			nicheX=$OPTARG;
			;;
		y)
			echo "The value of -y (spatial y) is $OPTARG";
			nicheY=$OPTARG;
			;;
		e)
			echo "The value of -e (ENM filename) is $OPTARG";
			nicheFile=$OPTARG;
			;; 
		c)
			echo "The value of -c (competition) is $OPTARG";
			competition=$OPTARG;
			;;
		i)
			echo "The value of -i (isolation by distance) is $OPTARG";
			ibd=$OPTARG;
			;;
		l)
			echo "The value of -l (use spatial) is $OPTARG";
			spatial=$OPTARG;
			;;
		t)
			echo "The Value of -t (treeSeq) is $OPTARG";
			treeSeq=$OPTARG;
			;;
		B)
			echo "The value of -B (number of subpops) is $OPTARG";
			subpopCount=$OPTARG;
			;;
		a) 
			echo "The value of -a (perform grow/select on all pops) is $OPTARG";
			allPops=$OPTARG;
			;;
			
		G) 
			echo "The value of -G (growth start gen) is $OPTARG";
			growthGen=$OPTARG;
			;;
		M) 
			echo "The value of -M (migration rate) is $OPTARG";
			migRate=$OPTARG;
			;;
		O) 
			echo "The value of -O (output gen) is $OPTARG";
			outputGen=$OPTARG;
			;;
		v) 
			echo "The value of -v (secondary contact yes/no) is $OPTARG";
			secContact=$OPTARG;
			;;
		V) 
			echo "The value of -V (secondary contact gen) is $OPTARG";
			secConGen=$OPTARG;
			;;
		H) 
			echo "The value of -H (selective sweep gen) is $OPTARG";
			sweepGen=$OPTARG;
			;;	
		X) 
			echo "The value of -X (scaling) is $OPTARG";
			scaling=$OPTARG;
			;;	
	esac
done

## generate the header using sample number and N
## where header is: "outSamp simulations"
header="${outSamp} ${simulations}"

echo "Begin running ${script}"
echo "Output to ${filename}"
echo "SLiM will output ${outSamp} samples"

echo
echo "~~~~~~~"
echo "HEADER:"
echo $header
echo "~~~~~~~"
echo

echo $header > $filename.header
echo " " >> $filename.header
#cat $filename.header

echo
echo "Beginning ${simulations} runs"
echo

for ((i=1;i<=simulations;i++)); do ## for run in simulations runs

	echo
	echo "#####"
	echo $i
	echo "#####"
	echo
	
	## for each iteration of slim,
	echo "Running SLiM"

	## note: you can pass extra constants to SLiM even if they aren't called by SLIM
	
	command="$slimPath -d N_un=$N -d mu_un=$mu -d recomb_un=$recomb \
	-d g_rate=$g_rate -d sweepType=\'$sweepType\' -d selStrength=$selStrength \
	-d outSamp=$outSamp -d N_dif_un=$N_dif -d nicheX=$nicheX \
	-d nicheY=$nicheY -d nicheFile=\'$nicheFile\' -d competition=$competition -d ibd=$ibd \
	-d spatial=$spatial -d treeSeq=$treeSeq -d subpopCount=$subpopCount \
	-d allPops=$allPops -d growthGen_un=$growthGen -d migRate=$migRate \
	-d outputGen_un=$outputGen -d secConGen_un=$secConGen -d secContact=$secContact \
	-d sweepGen_un=$sweepGen \
	-d scaling=$scaling \
	-d textName=\'$filename-$TIMESTAMP.txt\' \
	-d tempName=\'$filename-$TIMESTAMP-$i.temp\' \
	-d treeName=\'$filename-$TIMESTAMP-$i.trees\' \
	-d locsName=\'$filename-$TIMESTAMP-$i.locs\' \
	$path$script >> slim-$filename-$script-$TIMESTAMP.log 2>&1"

	echo "----------------------------"
	echo "Command is:"
	echo $command
	echo "----------------------------"

	## run the simulation and let it output as needed
	## will output all the runs to .txt and the last run to .temp (inside the slim script being called)
	
	eval $command

done 