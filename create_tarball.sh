#!/bin/bash

#define version number
version=10.3beta 

# define tarball name
tarball_tmp=flexpart_v$version 

# clean old package
rm -r $tarball_tmp

# create basic dir structure
mkdir $tarball_tmp
mkdir $tarball_tmp/src

# copy fortran source files
cp src/*.f90 $tarball_tmp/src
cp -r src/gributils $tarball_tmp/src 

# copy makefile
cp src/makefile $tarball_tmp/src
#cp src/makefile.gfs $tarball_tmp/src

# copy default options (for the distribution they work with the defult flex_ecmwf test winds)
# cp -r options $tarball_tmp
cp -r options_flex_ecmwf_EA $tarball_tmp/options

# copy default pathnames 
cp pathnames_distribution $tarball_tmp/pathnames


# add ECMWF retrieve routines
mkdir $tarball_tmp/preprocess
#mkdir $tarball_tmp/preprocess/flex_ecmwf
mkdir $tarball_tmp/preprocess/flex_extract
#cp -r flex_ecmwf_src/* $tarball_tmp/preprocess/flex_ecmwf/
cp -r flex_ecmwf_src/* $tarball_tmp/preprocess/flex_extract/

# copy default AVAILABLE
cp AVAILABLE_flex_ecmwf_EA $tarball_tmp/AVAILABLE


# directory for reading routines
mkdir $tarball_tmp/postprocess

# add fortran reading routines
mkdir $tarball_tmp/postprocess/flex_read_fortran
cp postprocess/flex_read_fortran/*.f $tarball_tmp/postprocess/flex_read_fortran
cp postprocess/flex_read_fortran/*.f90 $tarball_tmp/postprocess/flex_read_fortran
cp postprocess/flex_read_fortran/makefile $tarball_tmp/postprocess/flex_read_fortran

# add matlab reading routines
mkdir $tarball_tmp/postprocess/flex_read_matlab
cp postprocess/flex_read_matlab/*.m $tarball_tmp/postprocess/flex_read_matlab

# examples
cp -r examples $tarball_tmp/ 


mkdir $tarball_tmp/tests

cp -r tests/NILU/test_1 $tarball_tmp/tests/
cp -r tests/flex_gen_cases $tarball_tmp/tests/


#return

#tar -cvf flexpart$version.tar $tarball_tmp/*
echo now can run:  "tar -cvf $tarball_tmp.tar $tarball_tmp/*"
echo preliminary: scp flexpart_v10.3beta.tar njord:public_html/FLEXPART
