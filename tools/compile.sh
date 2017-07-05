#!/bin/bash

root="/snfs01/ou/models/pcom1.1"
platform="bcm.ifc"

source $root/bin/environs.$platform
exp_dir=$PWD
executable="$exp_dir/pcom.x"
#code_dir="$root/src"
code_dir="$exp_dir/my_src"

build_dir="$exp_dir/build"
mkmfTemplate="$root/bin/mkmf.template.$platform"
mkmf_exec="$root/bin/mkmf -f -m Makefile -a $code_dir -t $mkmfTemplate -p
$executable"
src_list=("pcom")

mkdir -p $build_dir
cd $build_dir

#make clean  # this line should be remove for regular compilation

$mkmf_exec $src_list

make

if [ $? -ne 0 ]; then
    echo "Make failed to create the executable $executable"
    exit 1
 else
    echo "Compilation done."
fi
