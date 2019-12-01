#!/bin/bash
# Ignacio Pisso, May 2017 

echo CREATE A NEW FLEXPART DISTRIBUTION

# get current commit hash
githash=$(git rev-parse --short --verify HEAD)
echo githash $githash 
#define version number with hash
version=10.4_$githash
echo version $version  

# define tarball name
targetdir=../flexpart_distribution/
echo targetdir $targetdir

# name distribution version
distribution_name=flexpart_v$version

# name distribution temporary dir
tarball_tmp=${targetdir}flexpart_v$version
echo tarball_tmp $tarball_tmp

# name distribution tarball file
#tarball=${targetdir}flexpart_v$version.tar 
tarball=${tarball_tmp}.tar 
echo tarball $tarball

# if needed clean old package
if [ -d $tarball_tmp ]; then
  echo
  echo clean old tarball
  hora=$(date +"%Y-%m-%d_%H%M%S")
  tarball_tmp_bk=$tarball_tmp$tarball_tmp_$hora
  echo tarball_tmp=$tarball_tmp exists: move to tarball_tmp_bk=$tarball_tmp_bk #and exit  
  mkdir $tarball_tmp_bk 
  mv $tarball_tmp $tarball_tmp_bk/ 
  mv $tarball $tarball_tmp_bk/ 
  #exit 
  echo old files moved to tarball_tmp_bk=$tarball_tmp_bk 
  echo
fi

# start packing


## needs in addition to the git repo ANCILLARY git repos
# VERIFY THESE RESOURCES EXIST BEFORE PACKING DISTRIBUTION
#1 OH file OH_variables.bin || OH_variables=../flexin/OH_FIELDS/OH_variables.bin
#2 flex_extract || flex_extract=../flex_extract_v7.0.4/
#3 flex_read_fortran from ../flex_read_matlab/export_basic  TODO: add functions to ../flex_read_matlab/export/ 
#  flex_read_matlab_src=../flex_read_matlab/export_basic
#4 tests/examples ../flex_tests_examples/examples3/*
#5  
#6




# mkdir container
echo ---------------------------------------------------------
echo ')' create basis dir $tarball_tmp
mkdir $tarball_tmp
echo ---------------------------------------------------------

echo

# patnames
echo ---------------------------------------------------------
echo ')' copy pathnames 
#cp pathnames_distribution $tarball_tmp/pathnames
cp pathnames $tarball_tmp/pathnames
echo ---------------------------------------------------------

echo 

# fortran source files
echo ---------------------------------------------------------
echo ')' copy src/
mkdir $tarball_tmp/src
cp src/*.f90 $tarball_tmp/src
cp -r src/gributils $tarball_tmp/src 
# echo '3)' copy makefile
cp src/makefile $tarball_tmp/src
#cp src/makefile.gfs $tarball_tmp/src
echo ---------------------------------------------------------

echo

# options dir
echo ---------------------------------------------------------
echo ')' copy options/ 
echo ---------------------------------------------------------
# (for the distribution they work with the defult flex_ecmwf test winds)
#cp -r options_flex_ecmwf_EA $tarball_tmp/options
mkdir $tarball_tmp/options
user_input_files="AGECLASSES     COMMAND        IGBP_int1.dat  OUTGRID        OUTGRID_NEST   RECEPTORS      RELEASES       surfdata.t     surfdepo.t"
for i in $user_input_files
do
  echo $i
  cp -r options/$i $tarball_tmp/options
  #echo copy $i to $tarball_tmp/options
done
mkdir $tarball_tmp/options/SPECIES
cp options/SPECIES/SPECIES* $tarball_tmp/options/SPECIES/
cp options/SPECIES/specoverview.f90 $tarball_tmp/options/SPECIES/
echo copy options/SPECIES/ to $tarball_tmp/options/SPECIES/
echo ---------------------------------------------------------

echo

# OH file 
echo ---------------------------------------------------------
echo ')' copy OH_variables.bin to flexin 
mkdir $tarball_tmp/flexin
OH_variables=../flexin/OH_FIELDS/OH_variables.bin
cp $OH_variables $tarball_tmp/flexin/
echo ---------------------------------------------------------

echo

# AVAILABLE
echo ---------------------------------------------------------
echo ')' copy AVAILABLE
#cp AVAILABLE_flex_ecmwf_EA $tarball_tmp/AVAILABLE
cp AVAILABLE $tarball_tmp/AVAILABLE
echo ---------------------------------------------------------

echo 

# output
echo ---------------------------------------------------------
echo  ')' create output/ #  mkdir $tarball_tmp/output
mkdir $tarball_tmp/output
echo ---------------------------------------------------------
echo output reference?
echo ---------------------------------------------------------

echo

# preprocess
echo ---------------------------------------------------------
echo ')' preprocess/
mkdir $tarball_tmp/preprocess
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
flex_extract=../flex_extract_v7.0.4/
echo include flex_extract v7.0.4 b7c1c04a204c91e53759ef590504bf52dfaece64
echo from $flex_extract [use git modules?] IP 3/2018 
cp $flex_extract/README.md $tarball_tmp/preprocess/flex_extract
cp -r $flex_extract/docs $tarball_tmp/preprocess/flex_extract
cp -r $flex_extract/grib_templates $tarball_tmp/preprocess/flex_extract
cp -r $flex_extract/python $tarball_tmp/preprocess/flex_extract
cp -r $flex_extract/src $tarball_tmp/preprocess/flex_extract
echo flex_extract copied
echo ---------------------------------------------------------
echo AVAILABLE generation scripts?
echo ---------------------------------------------------------

#echo '10)' cp example generating scripts [a separate repository]
#echo moved below
#mkdir $tarball_tmp/examples
#cp -r examples/*.sh $tarball_tmp/examples/ 
#cp -r examples/Makefile $tarball_tmp/examples/ 

echo 

# postprocess
echo ---------------------------------------------------------
echo ')' postprocess/

postprocess=postprocess
mkdir $tarball_tmp/$postprocess
echo -----------------flex_read_fortran-------------------
flex_read_fortran_src=$postprocess/flex_read_fortran/
#echo ')'  directory for reading routines
#echo '12)'  add fortran reading routines [a separate repository]
mkdir $tarball_tmp/$postprocess/flex_read_fortran
cp $postprocess/flex_read_fortran/*.f $tarball_tmp/$postprocess/flex_read_fortran
cp $postprocess/flex_read_fortran/*.f90 $tarball_tmp/$postprocess/flex_read_fortran
cp $postprocess/flex_read_fortran/makefile $tarball_tmp/$postprocess/flex_read_fortran
echo flex_read_fortran copied from $flex_read_fortran_src

echo -----------------flex_read_matlab-------------------
flex_read_matlab_src=../flex_read_matlab/export_basic
mkdir $tarball_tmp/$postprocess/flex_read_matlab
cp -r $flex_read_matlab_src/* $tarball_tmp/$postprocess/flex_read_matlab
echo flex_read_fortran from $flex_read_matlab_src  # NOT copied
# add matlab reading routines
#mkdir $tarball_tmp/postprocess/flex_read_matlab
#cp postprocess/flex_read_matlab/*.m $tarball_tmp/postprocess/flex_read_matlab
echo ---------------------------------------------------------

echo

echo ---------------------------------------------------------
echo ')' tests/
###############################################################
#echo '13) tests' 
mkdir $tarball_tmp/tests
#echo -----------------flex_read_fortran-------------------
#echo 'b) ./tests/flex_read_fortran/' 
#echo fixme
#mkdir $tarball_tmp/tests/flex_read_fortran
#cp tests/flex_read_fortran/test_read_default.sh  $tarball_tmp/tests/flex_read_fortran


###############################################################
echo ----------------- examples -------------------

#echo ') ./tests/examples/' 
mkdir $tarball_tmp/tests/examples
#echo ') scripts' 
#cp -r ./tests/examples/*.sh $tarball_tmp/tests/examples/ 
#echo ') makefile' 
#cp -r ./tests/examples/Makefile $tarball_tmp/tests/examples/ 
# echo USAGE: ~/repos/flexpart/tests/examples'$' make run
cp -r ../flex_tests_examples/examples3/* $tarball_tmp/tests/examples/

echo ----------------- examples_reference -------------------
cp -r ./tests/examples_reference $tarball_tmp/tests/



###############################################################
# echo -----------------postprocess examples-------------------
echo -----------------read examples-------------------

#echo '13 c) ./tests/read_examples/' 
mkdir $tarball_tmp/tests/read_examples

cp tests/read_examples/declare_examples $tarball_tmp/tests/read_examples/
cp tests/read_examples/display_examples.sh $tarball_tmp/tests/read_examples/
cp tests/read_examples/examples_output.txt $tarball_tmp/tests/read_examples/
cp tests/read_examples/read_examples.sh $tarball_tmp/tests/read_examples/
cp tests/read_examples/read_grids.sh $tarball_tmp/tests/read_examples/
cp tests/read_examples/read_parts.sh $tarball_tmp/tests/read_examples/
cp tests/read_examples/set_examples_all $tarball_tmp/tests/read_examples/
cp tests/read_examples/set_examples_3.sh $tarball_tmp/tests/read_examples/
cp tests/read_examples/read_headers.sh $tarball_tmp/tests/read_examples/
cp tests/read_examples/read_examples_output.txt $tarball_tmp/tests/ #read_examples/
# echo USAGE ~/repos/flexpart/tests/read_examples'$'./read_grids.sh

###############################################################
echo ------------compare examples-------------------
#echo tests/compare_examples.sh
mkdir $tarball_tmp/tests/compare_examples


#mkdir $tarball_tmp/tests/compare_examples
#cp tests/compare_examples/compare_grids.sh $tarball_tmp/tests/compare_examples
#cp tests/compare_grids.sh $tarball_tmp/tests/
cp tests/compare_examples/*.sh $tarball_tmp/tests/compare_examples
cp tests/compare_examples/compare_grids_output.txt $tarball_tmp/tests/ #compare_examples
# list of examples with units
#cp tests/declare_examples $tarball_tmp/tests/

echo 

# ~/repos/flexpart/tests$./compare_grids.sh 

#echo mkdir $tarball_tmp/tests/examples2/
#echo cp tests/examples2/setup.sh $tarball_tmp/tests/examples2/
# echo --repeat examples-------------------
# echo FIXME 

###############################################################
#echo -----------------ctbto-------------------
# mkdir $tarball_tmp/tests/ctbto

# cp -r tests/NILU/test_1 $tarball_tmp/tests/
# cp -r tests/default_cases $tarball_tmp/tests/

echo ---------------------------------------------------------
echo create tarball
#tar cvf $tarball  $tarball_tmp  
#tar cf $tarball  $tarball_tmp  
#cd 

cd $targetdir
tar cf $distribution_name.tar $distribution_name 

pwd


echo  tarball $tarball complete
echo exported untarred files in $tarball_tmp 

echo cp -r preprocess/flex_extract/work $tarball_tmp/preprocess/flex_extract/ 
echo cd $tarball_tmp/src
echo $HOME/repos/flexpart/src/make_in_laptop.sh
echo cd .. ';' ./src/FLEXPART 
echo cd postprocess/flex_read_fortran
echo make test
echo  max:  0.115784094     mean:   4.70877676E-05
#echo cd $tarball_tmp/tests/examples ';'   make run
echo cd ../../tests/examples ';'   make run
#echo cd $tarball_tmp/tests/read_examples
echo cd ../read_examples
echo ./read_examples.sh
echo ./read_examples.sh '>' ../read_examples_output.txt
echo cd ../compare_examples
echo ./compare_grids.sh
echo ./compare_grids.sh '>' ../compare_grids_output.txt

echo e.g. tar --append --file=$tarball_tmp/ ../compare_grids_output.txt ../read_examples_output.txt
 
 


exit
#return
###############################################################

# obtain $FLEXHOME (and set)
#1 cd $FLEXHOME/src 

#2 compile
#
#[laptop] source /Users/ignacio/repos/flexpart/src/make_in_laptop.sh 
# [njord] make
# ->created executable (FLEXPART)

#3 execute in src (absolute paths)
#
#[laptop] cp  /Users/ignacio/repos/flexpart/src/pathnames .
#[njord] FIXME
#
# mkdir output
# ./FLEXPART
# ->created output in output/

#4 read output
# cd  $FLEXHOME/postprocess/flex_read_fortran/
# make
# -> printheader* printgrid* flex_read_compare2*
#/postprocess/flex_read_fortran$./printheader ../../src/output/
#/postprocess/flex_read_fortran$./printgrid ../../src/output/ conc
# -> output in stdout (max:   11122924.0     sum:   90330784.0)

#5 execute in $FLEXHOME
# cd $FLEXHOME
# get winds
#[laptop] cp -r ~/repos/flex_winds/work/ ./preprocess/flex_extract/ 
#[njord] curl https://folk.nilu.no/~ignacio/FLEXPART/EA120101.tar --output EA120101.tar ; tar -xvf EA120101.tar ; mv flex_extract/work preprocess/flex_extract/ ; rmdir flex_extract

# src/FLEXPART
# -> output in $FLEXHOME/output/

#6 read output
# postprocess/flex_read_fortran/printheader output/
# postprocess/flex_read_fortran/printgrid output/ conc
# -> output in stdout ( max:   11578738.0     sum:   104058720.)

#7 gnererate examples
# cd $FLEXHOME/tests/examples

#make run

#make examples
#make batch
#./run_batch_cl.sh

#make (set_default_example.sh)
#tests/examples$../../src/FLEXPART
#output

#8 read examples:
#cd $FLEXHOME/tests/read_examples  
# ./read_headers.sh
# ./read_grids.sh 

#9 compare examples with reference
