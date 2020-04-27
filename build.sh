#!/usr/bin/env bash

builddir=build
srcdir=$(dirname $0)
n_jobs=1

function error() {
    echo
    echo "Got error< $1 >"
    echo
    exit 1
}

echo
echo Updating submodules...
echo

git submodule init
git submodule update

[ -d $builddir ] && rm -rf $builddir
mkdir -p $builddir
cmake_args=" \
    -DEigen3_DIR=$srcdir/eigen3 \
    "

pushd $builddir
cmake $cmake_args $srcdir || error 'CMake error'
make -j $n_jobs || error 'Error building project'
popd

echo
echo Successfully built project
echo
