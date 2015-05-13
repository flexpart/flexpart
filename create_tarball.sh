#!/bin/bash

version=9.2.0.2 
tarball_tmp=flexpart$version 

mkdir $tarball_tmp
mkdir $tarball_tmp/src
cp src/*.f90 $tarball_tmp/src
cp src/makefile $tarball_tmp/src
cp src/makefile.gfs $tarball_tmp/src
cp -r options $tarball_tmp

mkdir $tarball_tmp/tests
cp -r tests/NILU/test_1 $tarball_tmp/tests/

#return

tar -cvf flexpart$version.tar $tarball_tmp/*

