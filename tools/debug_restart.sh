#!/usr/bin/env python

import os
import sys
import glob

def runCmd(cmd):
  print(cmd)
  stat = os.system(cmd)
  if stat != 0:
    print("Error happen when run: "+cmd)
    sys.exit()

def runCmd_notstop(cmd):
  print(cmd)
  stat = os.system(cmd)

# convert fix time to record time
mydir = '/Users/ou/pcom/'

runCmd("make compile")
runCmd("cp namelist_norst namelist")
runCmd("make run")
runCmd("cp namelist_rst namelist")
runCmd("make run")

rsts   = glob.glob('output/debug/*_rst.*')
norsts = glob.glob('output/debug/*_norst.*')

for i in range(0,len(rsts)):
    cmd = 'cmp '+norsts[i]+' '+rsts[i]
    runCmd_notstop(cmd)
