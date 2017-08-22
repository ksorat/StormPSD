#! /bin/bash

#Moves unfinished blocks out of the way
OUTDIR="unfin"
INJSIZE="900M"
TRAPSIZE="100M"

COMS="ls -lh"
find Trap -name "*.h5part" -size -${TRAPSIZE} -exec ${COMS} {} +
find Inj* -name "*.h5part" -size -${INJSIZE}  -exec ${COMS} {} +