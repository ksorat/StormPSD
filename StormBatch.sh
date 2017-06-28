#!/bin/bash

#Uses arguments to submit batch of jobs
#StormBatch.sh B0 Bend

#Will submit jobs w/ runIDs in [B0,Bend]
#B0/Bend default to 1/10, "StormInj"

export B0=${1:-1}
export BEND=${2:-10}
export STUB=$3


export WALL="12:00" #Wallclock time
export NUMT="16" #Run with NUMT threads on 1 node

export QUEUE="regular"
export LOG="$STUB.%I.log" #Set names for log files

bsub -a poe -P "UJHB0003" -W $WALL -n 1 -q $QUEUE -J "$STUB[$B0-$BEND]" -e ${LOG} -o ${LOG} "StormBlock.sh"
