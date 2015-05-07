#!/bin/bash


tarball_tmp=flexpart$1 

mkdir $tarball_tmp
mkdir $tarball_tmp/src
cp src/*.f90 $tarball_tmp/src
cp src/makefile $tarball_tmp/src
cp src/makefile.gfs $tarball_tmp/src
cp -r options $tarball_tmp

tar -cvf flexpart$1.tar $tarball_tmp/*

