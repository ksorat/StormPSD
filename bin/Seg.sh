#! /bin/bash

#Moves unfinished blocks out of the way
INJSIZE="900M"
TRAPSIZE="100M"

#COMS="ls -lh"
COMS="mv -t unfin"
find Trap -name "*.h5part" -size -${TRAPSIZE} -exec ${COMS} {} +
find Inj* -name "*.h5part" -size -${INJSIZE}  -exec ${COMS} {} +