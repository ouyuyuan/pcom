#!/bin/bash

# debug by compare output with stand outputs

cd check_output/

standfiles=`ls standard_*.nc`

for file in $standfiles
do
  stand=$file
  debug="${stand/standard_/}"
  cmd="cmp $debug $stand"

  echo $cmd
  $cmd

#  if [ $? -ne 0 ]; then
#    exit 1
#  fi
done
