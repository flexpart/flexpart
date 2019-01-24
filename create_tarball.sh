#!/bin/bash

#define version number
version=10.3beta5

# define tarball name
targetdir=../flexpart_distribution/
tarball_tmp=${targetdir}flexpart_v$version 

# clean old package
if [ -d $tarball_tmp ]; then
  echo $tarball_tmp exists: move to $tarball_tmp.bk   
  mv $tarball_tmp ${tarball_tmp}.bk 
  exit 
fi

echo ---------------------------------------------------------
echo ')' create basic dir structure
mkdir $tarball_tmp
echo ---------------------------------------------------------
echo ---------------------------------------------------------
##############################################################
echo ')' pathnames 
#cp pathnames_distribution $tarball_tmp/pathnames
cp pathnames $tarball_tmp/pathnames
echo ---------------------------------------------------------
##############################################################
echo ')' src/
mkdir $tarball_tmp/src
cp src/*.f90 $tarball_tmp/src
cp -r src/gributils $tarball_tmp/src 
# echo '3)' copy makefile
cp src/makefile $tarball_tmp/src
#cp src/makefile.gfs $tarball_tmp/src
echo ---------------------------------------------------------
################################################################
echo ')' options 
# (for the distribution they work with the defult flex_ecmwf test winds)
#cp -r options_flex_ecmwf_EA $tarball_tmp/options
mkdir $tarball_tmp/options

user_input_files="AGECLASSES     COMMAND        IGBP_int1.dat  OUTGRID        OUTGRID_NEST   RECEPTORS      RELEASES       surfdata.t     surfdepo.t"

for i in $user_input_files
do
  echo $i
  cp -r options/$i $tarball_tmp/options
done


mkdir $tarball_tmp/options/SPECIES
cp options/SPECIES/SPECIES* $tarball_tmp/options/SPECIES/
cp options/SPECIES/specoverview.f90 $tarball_tmp/options/SPECIES/
echo ---------------------------------------------------------
################################################################
echo ')' AVAILABLE
#cp AVAILABLE_flex_ecmwf_EA $tarball_tmp/AVAILABLE
cp AVAILABLE $tarball_tmp/AVAILABLE

echo ---------------------------------------------------------
################################################################
echo  ')' output / #  mkdir $tarball_tmp/output
mkdir $tarball_tmp/output
echo ---------------------------------------------------------
################################################################
echo ')' preprocess/
mkdir $tarball_tmp/preprocess
#############################
echo -----------------flex_extract-------------------
#echo '6)'  mkdir $tarball_tmp/flex_extract [a separate repository]
#mkdir $tarball_tmp/preprocess
#mkdir $tarball_tmp/preprocess/flex_ecmwf
mkdir $tarball_tmp/preprocess/flex_extract

#echo '7)  add ECMWF retrieve routines (change EA wind files for latest source code)'
#mkdir $tarball_tmp/preprocess/flex_extract
#mkdir $tarball_tmp/preprocess/flex_extract/work
#cp -r flex_ecmwf_src/* $tarball_tmp/preprocess/flex_ecmwf/
#cp -r flex_ecmwf_src/* $tarball_tmp/preprocess/flex_extract/
## cp -r flex_extract/work/EA* $tarball_tmp/preprocess/flex_extract/work   




#echo '10)' cp example generating scripts [a separate repository]
#echo moved below

#mkdir $tarball_tmp/examples
#cp -r examples/*.sh $tarball_tmp/examples/ 
#cp -r examples/Makefile $tarball_tmp/examples/ 
echo ---------------------------------------------------------
################################################################
echo postprocess/

postprocess=postprocess
mkdir $tarball_tmp/$postprocess
echo -----------------flex_read_fortran-------------------
#echo ')'  directory for reading routines
#echo '12)'  add fortran reading routines [a separate repository]
mkdir $tarball_tmp/$postprocess/flex_read_fortran
cp $postprocess/flex_read_fortran/*.f $tarball_tmp/$postprocess/flex_read_fortran
cp $postprocess/flex_read_fortran/*.f90 $tarball_tmp/$postprocess/flex_read_fortran
cp $postprocess/flex_read_fortran/makefile $tarball_tmp/$postprocess/flex_read_fortran

echo -----------------flex_read_matlab-------------------

# add matlab reading routines
#mkdir $tarball_tmp/postprocess/flex_read_matlab
#cp postprocess/flex_read_matlab/*.m $tarball_tmp/postprocess/flex_read_matlab

###############################################################

echo ---------------------------------------------------------
echo tests/

#echo '13) tests' 
mkdir $tarball_tmp/tests

###############################################################
echo -----------------flex_read_fortran-------------------

#echo 'b) ./tests/flex_read_fortran/' 
echo fixme
#mkdir $tarball_tmp/tests/flex_read_fortran
#cp tests/flex_read_fortran/test_read_default.sh  $tarball_tmp/tests/flex_read_fortran

###############################################################
echo -----------------examples-------------------

#echo ') ./tests/examples/' 
mkdir $tarball_tmp/tests/examples
echo ') scripts' 
cp -r ./tests/examples/*.sh $tarball_tmp/tests/examples/ 
echo ') makefile' 
cp -r ./tests/examples/Makefile $tarball_tmp/tests/examples/ 

# echo USAGE: ~/repos/flexpart/tests/examples'$' make run


###############################################################
echo -----------------postprocess examples-------------------
echo --read examples-------------------

#echo '13 c) ./tests/read_examples/' 
mkdir $tarball_tmp/tests/read_examples
cp tests/read_examples/read_grids.sh $tarball_tmp/tests/read_examples/
cp tests/read_examples/read_headers.sh $tarball_tmp/tests/read_examples/

# echo USAGE ~/repos/flexpart/tests/read_examples'$'./read_grids.sh

###############################################################
echo --compare examples-------------------
#echo tests/compare_examples.sh

#mkdir $tarball_tmp/tests/compare_examples
#cp tests/compare_examples/compare_grids.sh $tarball_tmp/tests/compare_examples
cp tests/compare_grids.sh $tarball_tmp/tests/
# list of examples with units
cp tests/declare_examples $tarball_tmp/tests/



# ~/repos/flexpart/tests$./compare_grids.sh 

#echo mkdir $tarball_tmp/tests/examples2/
#echo cp tests/examples2/setup.sh $tarball_tmp/tests/examples2/
echo --repeat examples-------------------
#echo FIXME 

###############################################################
echo -----------------ctbto-------------------
mkdir $tarball_tmp/tests/ctbto

# cp -r tests/NILU/test_1 $tarball_tmp/tests/
# cp -r tests/default_cases $tarball_tmp/tests/


exit
#return
###############################################################



