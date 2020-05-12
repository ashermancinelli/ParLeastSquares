#!/bin/bash

builddir=build
[ -d $builddir ] && rm -rf $builddir
mkdir $builddir
cd $builddir

cmake .. -DEigen3_DIR=$INSTALL_LOCAL/share/eigen3/cmake || exit 1
make -j 8 || exit 1
ctest || exit 1

echo
echo Built and tested ParLeastSquares
echo
