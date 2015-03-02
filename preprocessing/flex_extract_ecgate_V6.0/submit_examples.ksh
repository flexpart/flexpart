#!/bin/ksh

echo 'submitting ECMWFDATA example scripts ..'

SERVER=`uname -n | cut -c 1-4`
if [[ ${SERVER} != ecgb ]] ; then
# submit via gateway from local server


for FILE in `ls CONTROL_ERA__*` ; do

  name=`echo $FILE | sed s,CONTROL_ERA__,,`
  echo submitting example  on demand script flex_ecmwf_$name
  if [ ! -f flex_ecmwf_$name ] ; then
    echo flex_ecmwf_$name does not exist
    echo run update_script.ksh before running this script.
    exit
  fi
  ecaccess-job-submit -queueName ecgb  flex_ecmwf_$name
done

echo submitting example  operational script ecmwf_idc_ops_ecgate
ecaccess-job-submit -queueName ecgb  ecmwf_idc_ops_ecgate

else
# submit on ecgate

for FILE in `ls CONTROL_ERA__*` ; do

  name=`echo $FILE | sed s,CONTROL_ERA__,,`
  echo submitting example  on demand script flex_ecmwf_$name
  if [ ! -f flex_ecmwf_$name ] ; then
    echo flex_ecmwf_$name does not exist
    echo run update_script.ksh before running this script.
    exit
  fi
  sbatch  flex_ecmwf_$name
done

echo submitting example  operational script flex_ecmwf_$name
sbatch  ecmwf_idc_ops_ecgate

fi
