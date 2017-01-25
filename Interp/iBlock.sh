#!/bin/bash
#
#Format, iBlock.sh <FSTUB> <MAP> <VTIDIR>

export FSTUB=$1
export MAP=$2
export VTIDIR=$3


module restore lfmtp
module list

export JOBID=$LSB_JOBINDEX
printf -v PID %"02d" $JOBID 


echo "Converting files ..."
echo "Stub = $FSTUB"
echo "Map = $MAP"
echo "VTI Directory = $VTIDIR"

ls -lh ${FSTUB}${PID}-??-??Z.hdf

lfm2cart.py $MAP ${FSTUB}${PID}-??-??Z.hdf
echo "Finished converting files ..."
echo "Moving files ..."
mv -v ${FSTUB}${PID}-??-??Z.vti $VTIDIR

