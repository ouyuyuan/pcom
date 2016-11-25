#!/bin/bash

# get the compiled code from 101

dir101="ou@172.16.0.101:/snfs01/ou/models/pcom1.1/"
dirloc="/home/ou/archive/pcom_src/v1.1/"

rsync -avz --exclude-from 'exclude.txt' $dir101 $dirloc
