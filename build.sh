#!/bin/bash

do_install=false

message()
{
  echo -----\\
  echo "------> $*"
  echo -----/
}

while [[ $# -gt 0 ]]
do
  case $1 in
    --install)
      do_install=true
      shift
      ;;
    *)
      echo Command not found.
      exit 1
      ;;
  esac
done

build()
{
  message Building ParLeastSquares

  builddir=build
  [ -d $builddir ] && rm -rf $builddir
  mkdir $builddir
  pushd $builddir

  cmake \
    -DCMAKE_BUILD_TYPE=Debug \
    -DEigen3_DIR=$INSTALL_LOCAL/share/eigen3/cmake \
    .. || exit 1
  make -j 8 || exit 1
  # ctest || exit 1
  popd

  echo
  echo Built and tested ParLeastSquares
  echo
}

install()
{
  message Installing python module
  python -m pip install --user -e .
}

build
if $do_install
then
  install
fi

message Build successful
