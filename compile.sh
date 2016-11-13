#!/bin/csh
#===========================================================
# run script for the p-sigma coordinate model
#===========================================================
#setenv PATH /bin:/usr/bin:/home/zhangyu/openmpizhy/bin:/home/zhangyu/pgizhy/linux86-64/9.0-1/bin 
set filepath   = ./FILE
set targer     = pcom.out
set COMPILER_OPTIONS =  "-O2 -r8"
/bin/rm $targer
cd $filepath
/bin/rm *.[o]
#/disk1/soft/intel/impi/3.1/bin64/mpiifort $COMPILER_OPTIONS -c *.f90
/disk5/home/lyh/openmpi/bin/mpif90 $COMPILER_OPTIONS -c *.f90
/disk5/home/lyh/openmpi/bin/mpif90 $COMPILER_OPTIONS *.o -o $targer
#/disk1/soft/intel/impi/3.1/bin64/mpiifort $COMPILER_OPTIONS *.o -o $targer
mv $targer ../
