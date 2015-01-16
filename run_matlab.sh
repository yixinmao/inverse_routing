#!/bin/bash

# basin="pend"

basin=$1
nst=$2
sst=$3
spkt=$4
# ./run_route_assim.sh /usr/local/MATLAB/R2013a $basin $nst $sst $spkt
./run_route_assim.sh /usr/local/MATLAB/R2014b/ $basin $nst $sst $spkt

