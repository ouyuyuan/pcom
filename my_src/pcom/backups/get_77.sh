#!/bin/bash

# get the compiled code from 77

dir77="lyh@172.16.0.77:/disk5/home/lyh/models/pcom1.1/"
dirloc="/home/ou/archive/pcom_src/v1.1/"

rsync -avz --exclude-from 'exclude.txt' $dir77 $dirloc
