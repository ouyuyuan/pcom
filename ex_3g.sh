#!/bin/csh
#$ -S /bin/csh
#$ -cwd
#$ -j y
#$ -N ex_3g
#$ -pe make 40
#$ -V

/snfs01/zhangyu/openmpi/bin/mpirun -n 40 /snfs01/zhangyu/pcom/ex_3g/pcom.out >& /snfs01/zhangyu/pcom/ex_3g/pcom.log
