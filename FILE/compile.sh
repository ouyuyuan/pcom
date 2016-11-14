#!/bin/csh
#===========================================================
# compile script for the p-sigma coordinate model
#===========================================================
setenv PATH /bin:/usr/bin:/home/zhangyu/lacie/zhangyu/openmpi/bin:home/zhangyu/lacie/zhangyu/pgi/linux86-64/9.0-1/bin 
set filepath   = ./FILE
set targer     = pcom.out
set COMPILER_OPTIONS = "-Kieee -fast -r8"
set NETCDF_INC = "~/lacie/zhangyu/netcdf/include/"
set NETCDF_LIB = "~/lacie/zhangyu/netcdf/lib/"
set FCOMPILER = ~/lacie/zhangyu/openmpi/bin/mpif90
/bin/rm $targer
cd $filepath
/bin/rm *.[o]
$FCOMPILER $COMPILER_OPTIONS -c *.f90 -I $NETCDF_INC -L $NETCDF_LIB -lnetcdf
$FCOMPILER $COMPILER_OPTIONS *.o -o $targer -I $NETCDF_INC -L $NETCDF_LIB -lnetcdf
mv $targer ../
