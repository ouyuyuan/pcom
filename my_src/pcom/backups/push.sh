#!/bin/bash

# push the craft code to 101
# remember to pull from 101 again after sussessful compile

dir101="ou@172.16.0.101:/snfs01/ou/models/pcom1.1/"
dirloc="/home/ou/archive/pcom_src/v1.1/"

rsync -avz --exclude-from 'exclude.txt' $dirloc $dir101
