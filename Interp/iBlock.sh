#!/bin/bash
#
#Format, iBlock.sh <FSTUB> <MAP> <LFMDIR> <VTIDIR>

export FSTUB=$1
export MAP=$2
export LFMDIR=$3
export VTIDIR=$4

module restore lfmtp
module list

#Subtract 1 to remap 1-24 to 0-23
export JOBID=$(echo $(($LSB_JOBINDEX-1)))

printf -v PID %"02d" $JOBID 
echo "Using ID $PID"

echo "Converting files ..."
echo "Stub = $FSTUB"
echo "Map = $MAP"
echo "LFM HDF Directory = $LFMDIR"

ls -lh ${LFMDIR}/${FSTUB}${PID}-??-??Z.hdf

lfm2cart.py $MAP ${LFMDIR}/${FSTUB}${PID}-??-??Z.hdf
echo "Finished converting files ..."

echo "Moving files ..."
mv -v ${LFMDIR}/${FSTUB}${PID}-??-??Z.vti $VTIDIR

