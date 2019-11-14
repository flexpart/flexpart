#!/bin/bash
# Ignacio Pisso, May 2017 
# Changes 2018-2019

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
#4 tests/examples ../flex_tests_examples/examples/*




# mkdir container
echo ---------------------------------------------------------
echo ')' create basis dir $tarball_tmp
mkdir $tarball_tmp
echo ---------------------------------------------------------

echo

# patnames
echo ---------------------------------------------------------
echo ')' copy pathnames 
cp pathnames $tarball_tmp/pathnames
echo ---------------------------------------------------------

echo 

# fortran source files
echo ---------------------------------------------------------
echo ')' copy src/
mkdir $tarball_tmp/src
cp src/*.f90 $tarball_tmp/src
cp -r src/gributils $tarball_tmp/src 
cp src/makefile $tarball_tmp/src
echo ---------------------------------------------------------

cp LICENSE $tarball_tmp/LICENSE_GPLv3
cp src/flexpart_license.txt  $tarball_tmp/src

echo

# options dir
echo ---------------------------------------------------------
echo ')' copy options/ 
echo ---------------------------------------------------------
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
mkdir $tarball_tmp/preprocess/flex_extract
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

echo 

# postprocess
echo ---------------------------------------------------------
echo ')' postprocess/

postprocess=postprocess
mkdir $tarball_tmp/$postprocess
echo -----------------flex_read_fortran-------------------
flex_read_fortran_src=$postprocess/flex_read_fortran/
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
echo ---------------------------------------------------------

echo

echo ---------------------------------------------------------
echo ')' tests/
#echo '13) tests' 
mkdir $tarball_tmp/tests
echo ----------------- examples -------------------

mkdir $tarball_tmp/tests/examples
cp -r ../flex_tests_examples/examples3/* $tarball_tmp/tests/examples/

echo ----------------- examples_reference -------------------
cp -r ./tests/examples_reference $tarball_tmp/tests/

echo -----------------read examples-------------------

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

echo ------------compare examples-------------------
mkdir $tarball_tmp/tests/compare_examples


cp tests/compare_examples/*.sh $tarball_tmp/tests/compare_examples
cp tests/compare_examples/compare_grids_output.txt $tarball_tmp/tests/ #compare_examples

echo 

echo ---------------------------------------------------------
echo create tarball

cd $targetdir
tar cf $distribution_name.tar $distribution_name 

pwd

echo  tarball $tarball complete
echo exported untarred files in $tarball_tmp 

echo verify:
echo cp -r preprocess/flex_extract/work $tarball_tmp/preprocess/flex_extract/ 
echo cd $tarball_tmp/src
echo $HOME/repos/flexpart/src/make_in_laptop.sh
echo cd .. ';' ./src/FLEXPART 
echo cd postprocess/flex_read_fortran
echo make test
echo e.g.: max:  0.115784094     mean:   4.70877676E-05
echo cd ../../tests/examples ';'   make run
echo cd ../read_examples
echo ./read_examples.sh
echo ./read_examples.sh '>' ../read_examples_output.txt
echo cd ../compare_examples
echo ./compare_grids.sh
echo ./compare_grids.sh '>' ../compare_grids_output.txt

echo e.g. tar --append --file=$tarball_tmp/ ../compare_grids_output.txt ../read_examples_output.txt
 
