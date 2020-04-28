#!/bin/bash

source /etc/profile.d/modules.sh

builddir=build
srcdir=$(realpath $(dirname $0))
n_jobs=1
python=3.7.0
gcc=9.1.0
cmake=3.15.3
cuda=10.1.243

module purge
module load python/$python
module load gcc/$gcc
module load cmake/$cmake
module load cuda/$cuda

echo
echo Updating submodules...
echo

git submodule init
git submodule update

[ -d $builddir ] && rm -rf $builddir
mkdir -p $builddir

cmake_args=""

pushd $builddir
cmake $cmake_args $srcdir || error 'CMake error'
make -j $n_jobs || error 'Error building project'
popd

echo
echo Successfully built project
echo
