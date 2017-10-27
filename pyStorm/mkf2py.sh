#!/bin/bash

#Compiles using f2py
FOPT="-O3"
FFLAGS="-qopenmp"
FILE="kSmooth"
COMP="intelem"
#f2py --no-lower --opt="${FOPT}" --f90flags="${FFLAGS}" -lgomp -c -m $FILE ${FILE}.F90
f2py --fcompiler=${COMP} --no-lower --opt="${FOPT}" --f90flags="${FFLAGS}" -liomp5 -c -m $FILE ${FILE}.F90

