#!/bin/bash
#Runs one storm block

module restore lfmtp
module list

export OMP_NUM_THREADS=16
export JOBID=${LSB_JOBINDEX}

printf -v PID %"04d" $JOBID

export INPDECK="Inps/${STUB}.$PID.xml"

unset MP_PE_AFFINITY
export MP_TASK_AFFINITY=core:$OMP_NUM_THREADS
export KMP_STACKSIZE=64M
export OMP_STACKSIZE=64M

echo "Running w/ input deck $INPDECK and RunID $JOBID"
echo "Using $OMP_NUM_THREADS OMP Threads w/ affinity $MP_TASK_AFFINITY"
echo "Running on host `hostname` on `date`"
echo "Running command ..."
echo "<push.x $INPDECK $JOBID>"
push.x $INPDECK $JOBID
echo "Finished run on `date`"