#!/bin/bash
source /disk5/home/lyh/.bashrc
source /disk5/home/lyh/.bash_profile
cd /disk5/home/lyh/pcom/
export I_MPI_PIN=off
#/disk1/soft/intel/impi/3.2.0.011/bin64/mpirun -r ssh -n 8 /disk5/home/lyh/pcom/pcom.out
/disk5/home/lyh/openmpi/bin/mpirun -n 8 /disk5/home/lyh/pcom/pcom.out
