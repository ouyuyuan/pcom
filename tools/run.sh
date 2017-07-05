#!/bin/bash

# Calc CPUs needed from namelist
line=`cat namelist | grep -w npx`
npx=`echo $line | sed 's/[^0-9]*//g'`
line=`cat namelist | grep -w npy`
npy=`echo $line | sed 's/[^0-9]*//g'`

np=$(($npx*$npy))

while :
do
   result=`select_hosts $np`
   if [ $result == "failed" ]; then
      echo "Not enough available hosts for $np CPUs. Will try again after 5 minute..."
      sleep 300
   elif [ $result == "done" ]; then
      break
   else
      echo "Unprocessed output of select_hosts: $result"
      exit 1
   fi
done

mpirun -r ssh -n $np ./pcom.x
