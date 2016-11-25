#!/bin/bash

# push the craft code to 77
# remember to pull from 77 again after sussessful compile

dir77="lyh@172.16.0.77:/disk5/home/lyh/models/pcom1.1/"
dirloc="/home/ou/archive/pcom_src/v1.1/"

rsync -avz --exclude-from 'exclude.txt' $dirloc $dir77
