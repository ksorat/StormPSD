#!/bin/bash

#Uses arguments to submit batch of jobs
#Format, lfmtpBatch.sh <STUB> B0 Bend
#Format, iBatch.sh D0 DEnd

export NSTUB="StormI"

export TAG="LR60-Quad-15s-AEH_mhd"
export DATESTUB="2013-03"
export MAP="$HOME/Work/maps/Quad2Hirez.h5"
export LFMDIR="$HOME/Work/StormPSD/StormPSD_Data/lfmData/"
export VTIDIR="$HOME/Work/StormPSD/StormPSD_Data/"
export D0=${1:-0}
export DEND=${2:-0}

#Job parameters
export WALL="8:00" #Wallclock time
export QUEUE="regular"

for d in $(seq $D0 $DEND)
do

	printf -v DD %"02d" $d
	export FSTUB=${TAG}_$DATESTUB-${DD}T
	export LOG="Interp.$DD.%I.log" #Set names for log files
	export NAME="${NSTUB}${DD}"
	echo "Submitting $NAME to work on $FSTUB"
	bsub -a poe -P "UJHB0003" -W $WALL -n 1 -q $QUEUE -J "${NAME}[1-24]" -e ${LOG} -o ${LOG} "iBlock.sh $FSTUB $MAP $LFMDIR $VTIDIR"
done
