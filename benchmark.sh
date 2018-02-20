#!/bin/bash

setup_subdir ()
{
	if [ -d $1 ]; then
		echo "Subdirectory "$1" already exists!"
		return
	else
		mkdir $1
		#cd $1
		ln -s ~/Code/Spike/defaults.m $1/.
		ln -s ~/Code/Spike/benchmark.m $1/.
		ln -s $DIR/$EXEC $1/.
		cp ~/bin/SpikeNet/wdir/random_seeds.rsd $1/.
		#cd ..
		echo "Subdirectory "$1" set up!"
	fi
}

#DIR=/SpikeNet/build/Release
DIR=~/bin/SpikeNet/Debug
EXEC="Spike"
ARGS=" -rf benchmark.m"
#SEED=123
THREADS=$1
DYNAMIC="FALSE" # Set to true for xgrid simulations
SCHEDULING="static" 	# dynamic[, n]
			# guided[, n]
			# runtime
			# static[, n] (Default)
#HASHS="E0L0ExcitSpikes.dat E0L1ExcitSpikes.dat L1weightsEfE.dat L1ExcitSpikes.dat L0InhibSpikes.dat L1InhibSpikes.dat"
#L1affNeuronsEfE.dat ptL1weightsEfE.dat ptL1ExcitSpikes.dat

export OMP_NUM_THREADS=$THREADS
#printf "Set number of threads to %d" $1
export OMP_DYNAMIC=$DYNAMIC
export OMP_SCHEDULE=$SCHEDULING
#export GSL_RNG_SEED=$SEED

PARENT=`pwd`
rm -R 0
rm -R 1

setup_subdir 0
cd 0
rm *.dat *.tbz

# Run with benchmark parameters
$DIR/$EXEC$ARGS

#printf "\n"
#for hash in $HASHS
#do
#    md5 $hash
#done
#printf "\n"

hashOrig=`md5 datHashs`

cd $PARENT
if [ $# -gt 1 ]; then
	#cd $PARENT
	setup_subdir 1
	cd 1
	rm *.dat *.tbz

	NEWARGS=""
	shift
	# Run with new arguments
	while (( "$#" )); do
		NEWARGS+=" -p ${1}"
		shift
	done
	# Set of cli args: $@

	# Rerun with new arguments
	$DIR/$EXEC$ARGS$NEWARGS

	hashNew=`md5 datHashs`

	cd $PARENT

	if [ "$hashOrig" != "$hashNew" ]; then
		echo "MD5 differences recorded in hashDiffs"
		diff 0/datHashs 1/datHashs > hashDiffs
	fi

#else
#	cd $PARENT
fi

# Call Matlab for analysis : http://www.mathworks.com/help/techdoc/matlab_env/f8-4994.html;jsessionid=bBtNNwnGHwKlZ28zJpDvLcyzSHXVHLsqZJH3d1W2n10nVYsvvvn7!1423304219
#Matlab -nosplash -nodesktop -r "SpikeAnalysis(\'$PARENT/\'); exit;"
#stty echo
