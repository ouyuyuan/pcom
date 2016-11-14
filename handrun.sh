#!/bin/bash

# Calc CPUs needed from namelist
line=`cat namelist | grep -w npx`
npx=`echo $line | sed 's/[^0-9]*//g'`
line=`cat namelist | grep -w npy`
npy=`echo $line | sed 's/[^0-9]*//g'`

np=$(($npx*$npy))

#mpirun -r ssh -n $np ./pcom.x
mpirun -n $np ./pcom.x
