#!/bin/bash
rm -fr build
mkdir build
cd build
cmake -DBOOST_ROOT=/usr/include/boost148 -DBOOST_LIBRARYDIR=/usr/lib64/boost148 ../
make
$HOME/TASC/build/TASC -b -v -y $HOME/TASC/example/y.txt \
-x $HOME/TASC/example/x.txt \
-k $HOME/TASC/example/ercc.txt \
-o $HOME/TASC/example/out
