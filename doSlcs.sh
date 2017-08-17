#!/bin/bash

export EB="ebtab.txt"
export H5P="Data/Trap/StormTrap.0001.h5part"

fstSlcs.py --ns 0   --ne 100   $EB $H5P &
fstSlcs.py --ns 101 --ne 200   $EB $H5P &
fstSlcs.py --ns 201 --ne 300   $EB $H5P &
fstSlcs.py --ns 301 --ne 400   $EB $H5P &
fstSlcs.py --ns 401 --ne 500   $EB $H5P &
fstSlcs.py --ns 501 --ne 600   $EB $H5P &
fstSlcs.py --ns 601 --ne 700   $EB $H5P &
fstSlcs.py --ns 701 --ne 800   $EB $H5P &
fstSlcs.py --ns 801 --ne 900   $EB $H5P &
fstSlcs.py --ns 901 --ne 1000  $EB $H5P &
fstSlcs.py --ns 1001 --ne 1200 $EB $H5P &

