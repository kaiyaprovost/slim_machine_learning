#!/bin/bash

TIMESTAMP=`date +"%s"`

echo "TIMESTAMP:" 
echo $TIMESTAMP

# bash run_slim_sims_multiple_times.sh -S 100 -I demography_only.slim -O neutral_demog >> run_slim_neutral.log 2>&1 &
# bash run_slim_sims_multiple_times.sh -S 100 -I demography_only.slim -O growth_demog -g 1.03 -N 1000 -f 10000 >> run_slim_growth.log 2>&1 &

## set default arguments here 
#simulations=2 # number of reps to run 
#script="model.slim" # which slim file 
#filename="model" # used to specify outfile. ## something like p1_ms_decline for .txt and .temp 
##	if you have multiple filenames they can be in quotes but you must add a for loop
## this really only matters if you're running multiple populations though 
## later may need a Npop argument

#mu="1e-7" ## mutation rate 
#THETA=-1 ## shouldn't be set unless manually (in case you want to set theta instead of mu)
## temporary ignoring THETA
#recomb="1e-8" #recomb rate # default 1e-8
#N=1000 #beginning pop # default 1000
#N_dif=$N
#g_rate=1 #growth or contract # default 1 (no growth or contract)
#outSamp=20 # number genomes to output # default 20
#slimPath=slim
#selStrength=0
#sweepType="n"

## defining arguments!!!

# -a flag # use all pops or not 		# default 0 
# -B flag # number of subpops 			# default 1 
# -c flag # competition yes/no			# default 1 -- only works if spatial 1
# -e flag # ENM filename 				# default /Users/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/vireotest_slim_27x45.asc
# -f flag # ending pop size 			# default same as $N
# -G flag # growth generation start		# default 1000
# -g flag # growth/contract rate 		# default 1 (no growth or contract), under 1 shrinks, over 1 grows
# -H flag # sweep generation start		# default 1000
# -h flag # what type of sweep wanted 	# default "n" no sweep, "h" hard sweep, "s" soft sweep ## consider adding things for balancing, etc
# -i flag # ibd yes/no					# default 1 -- only works if spatial 1
# -I flag # which slim file 			# default "model.slim"
# -l flag # spatial yes/no 				# default 0
# -M flag # migration rate 				# default 0.1
# -m flag # mutation rate 				# default 1e-7
# -N flag # beginning pop size 			# default 1000
# -o flag # name of output file 		# default "modeltest"
# -O flag # output generation 			# default 2000
# -P flag # path 						# default 1 
# -p flag # where slim is located 		# default "slim" 
# -r flag # recombination rate 			# default 1e-8
# -s flag # number genomes to output 	# default 20
# -S flag # number of reps to run 		# default 2
# -t flag # use treeseq or not 			# default 1
# -v flag # secondary contact or not 	# default 0
# -V flag # secondary contact start		# default 1000
# -w flag # selection coefficient		# default 0.5 (beneficial), 0 neutral, negative deleterious
# -x flag # niche model X dimension		# default 27
# -y flag # niche model Y dimension		# default 45


while getopts p:S:I:P:o:m:r:N:g:f:s:w:h:x:y:e:c:i:l:t:B:a:G:M:O:v:V:H: OPTION
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
	esac
done


if [ ! -d "slim-$filename-$TIMESTAMP/" ]; then
  mkdir "slim-$filename-$TIMESTAMP/"
fi
echo "moving to slim-$filename-$TIMESTAMP/"
cd "slim-$filename-$TIMESTAMP/"

## set some arguments after because needed

## do checks of arguments to make sure going to be okay
## problem is that bash can't handle floats 
## may need to figure this out inside slim or something

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
	#slim $script;

	## string together all of the commands with -d, the only thing you need to do is toggle if theta or not 
	## note: you can pass extra constants to SLiM even if they aren't called by SLIM
	
	## THIS ONLY WORKS WITH SINGLE FILENAMES

	#command="$slimPath -d N=$N -d mu=$mu -d recomb=$recomb -d g_rate=$g_rate -d sweepType=\'$sweepType\' -d selStrength=$selStrength -d outSamp=$outSamp -d N_dif=$N_dif -d textName=\'$filename.txt\' -d tempName=\'$filename-$TIMESTAMP.temp\' $script >> slim-$filename-$script-$TIMESTAMP.log 2>&1"
	command="$slimPath -d N=$N -d mu=$mu -d recomb=$recomb -d g_rate=$g_rate -d sweepType=\'$sweepType\' -d selStrength=$selStrength -d outSamp=$outSamp -d N_dif=$N_dif -d nicheX=$nicheX -d nicheY=$nicheY -d nicheFile=\'$nicheFile\' -d competition=$competition -d ibd=$ibd -d spatial=$spatial -d treeSeq=$treeSeq -d subpopCount=$subpopCount -d allPops=$allPops -d growthGen=$growthGen -d migRate=$migRate -d outputGen=$outputGen -d secConGen=$secConGen -d secContact=$secContact -d sweepGen=$sweepGen -d textName=\'$filename-$TIMESTAMP.txt\' -d tempName=\'$filename-$TIMESTAMP-$i.temp\' -d treeName=\'$filename-$TIMESTAMP-$i.trees\' $path$script >> slim-$filename-$script-$TIMESTAMP.log 2>&1"

	echo "----------------------------"
	echo "Command is:"
	echo $command
	echo "----------------------------"

	## run the simulation and let it output as needed
	## will output all the runs to .txt and the last run to .temp (inside the slim script being called)
	
	eval $command
	
	if [ $treeSeq==1 ]; then
		echo "Generated $filename-$TIMESTAMP-$i.trees for simulation $i"
		
		echo "python -- adding mutations and outputting vcf file "
		source activate py36
		
		command2="python /Users/kprovost/Documents/Github/slim_machine_learning/Overlaying_neutral_mutations_klp.py $filename-$TIMESTAMP-$i.trees $mu $TIMESTAMP"
		
		echo "----------------------------"
		echo "Command2 is:"
		echo $command2
		echo "----------------------------"
		
		eval $command2
		
		echo "done" 
		
		echo "slim -- converting to MS format" 
		
		command3="perl /Users/kprovost/Documents/Github/slim_machine_learning/vcf2MS.pl $filename-$TIMESTAMP-$i-overlaid.vcf $filename-$TIMESTAMP-$i.temp $outSamp"
		
		echo "----------------------------"
		echo "Command3 is:"
		echo $command2
		echo "----------------------------"
		
		eval $command3
		
		echo "done"
	
		echo "proceeding to samplestats step" 
	
	fi
	
	
	echo "Generated $filename-$TIMESTAMP-$i.temp for simulation $i"
	echo "Adding to $filename-$TIMESTAMP.txt"

	if [ -e "$filename-$TIMESTAMP.txt" ]; then    
		cat $filename-$TIMESTAMP-$i.temp >> $filename-$TIMESTAMP.txt
	else
		cat $filename.header $filename-$TIMESTAMP-$i.temp > $filename-$TIMESTAMP.txt
	fi

	#ls *.temp

	#less $filename-$TIMESTAMP.txt

	## send to sample stats only if not doing treeSeq

	## cat the .temp file to sample_stat with the header on it 
	## append (>>) the sample_stat stuff to a new document
	echo "Sending to sample_stats"

	for f in $filename; do
		echo "(Output to ${f}-${TIMESTAMP}.sumstats.stats)"
		echo
		cat $filename.header $f-$TIMESTAMP-$i.temp > $f-$TIMESTAMP-$i.withheader.temp

		#cat $f-$TIMESTAMP.withheader.temp | /home/kprovost/nas1/Analysis_SLiM/msdir/sample_stats >> $f-$TIMESTAMP.sumstats.stats 
		
		########################################################
		## HERE IS COMMENTED OUT BECAUSE WON'T WORK ON HUXLEY ##
		########################################################
		
		#cat $f-$TIMESTAMP-$i.withheader.temp | /Users/kprovost/Documents/Classes/Machine_Learning/popGenMachineLearningExamples-master/msdir/sample_stats >> $f-$TIMESTAMP.sumstats.stats 
		#echo " " >> $f-$TIMESTAMP.sumstats.stats 

	done

done 

## after all simulations, call the python script 
#ls *-$TIMESTAMP.sumstats.stats | xargs -t python /home/kprovost/Documents/Classes/Machine_Learning/SLiMTreeSeqPub-master/models_to_run/slim_model_selection.py