#!/bin/ksh

if [[ $# -ne 4 && (-z "$GATEWAY" || -z "$DESTINATION" ||  -z "$ECUID" ||  -z "$ECGID" )  ]] ; then

  echo Syntax: update_script.ksh Gateway Destination EC_user_ID EC_group_ID
  echo e.g. update_script.ksh srvx7.img.univie.ac.at leo@genericSftp lh0 spatlh00
  echo or set environment variables  GATEWAY DESTINATION ECUID ECGID
exit
else
  if [[ $# -eq 4 ]] ; then
    export GATEWAY=$1
    export DESTINATION=$2
    export ECUID=$3
    export ECGID=$4
  fi
fi

########### Just two examples  ###############
# Please modify for your environment         #
##############################################

#[ -z "$GATEWAY" ] && GATEWAY=srvx7.img.univie.ac.at
#[ -z "$DESTINATION" ] && DESTINATION=leo@genericSftp
#[ -z "$ECUID" ] && ECUID=lh0
#[ -z "$ECGID" ] && ECGID=spatlh00

#[ -z "$GATEWAY" ] && GATEWAY=ctbto4.ctbto.org
#[ -z "$DESTINATION" ] && DESTINATION=atmops@ops
#[ -z "$ECUID" ] && ECUID=cbb

export VERSION=6
export SUBVERSION=0
echo 'ECMWFDATA V'${VERSION}.${SUBVERSION}' is installed for:'
echo Gateway     $GATEWAY 
echo Destination $DESTINATION 
echo ECMWF user ID $ECUID 
echo ECMWF group ID $ECGID 
echo Output will be written into '$SCRATCH' directory of user $ECUID - /scratch/ms/$ECGID/$ECUID

echo 'Note: These settings can be changed via environment variables'
echo '$GATEWAY $DESTINATION  $ECUID $ECGID'

cat flex_ecmwf_header_template | sed "s,xxx,${ECUID},"  | sed "s,ggg,${ECGID}," >flex_ecmwf_header

cat CONTROL_OPS_TEMPLATE | sed "s,xxx.xxx.xxx.xxx,${GATEWAY}," | sed "s,xxx@xxx,${DESTINATION},"  | sed "s,xxx,${ECUID},"  | sed "s,ggg,${ECGID},"  | sed "s,v.v,${VERSION}.${SUBVERSION}," > CONTROL_OPS_V${VERSION}.${SUBVERSION}
cat CONTROL_OPS_TEMPLATE | sed "s,xxx.xxx.xxx.xxx,${GATEWAY}," | sed "s,xxx@xxx,${DESTINATION},"  | sed "s,xxx,${ECUID},"  | sed "s,ggg,${ECGID},"  | sed "s,v.v,${VERSION}.${SUBVERSION}," | sed "s,M_TYPE AN FC FC FC FC FC FC FC FC FC FC FC AN FC FC FC FC FC FC FC FC FC FC FC 12,M_TYPE AN FC FC FC FC FC AN FC FC 4V FC FC 4V FC FC FC FC FC AN FC FC 4V FC FC 4V," | sed "s,M_TIME 00 00 00 00 00 00 00 00 00 00 00 00 12 12 12 12 12 12 12 12 12 12 12 12 12,M_TIME 00 00 00 00 00 00 06 00 00 09 00 00 09 12 12 12 12 12 18 12 12 21 12 12 21," | sed "s,M_STEP 00 01 02 03 04 05 06 07 08 09 10 11 00 01 02 03 04 05 06 07 08 09 10 11 12,M_STEP 00 01 02 03 04 05 00 07 08 00 10 11 03 01 02 03 04 05 00 07 08 00 10 11 03," | sed "s,DTIME 1,DTIME 3," | sed "s,M_ETA 1,M_ETA 0," | sed "s,M_GAUSS 0,M_GAUSS 1," | sed "s,M_SMOOTH 0,M_SMOOTH 179," > CONTROL_OPS_V${VERSION}.${SUBVERSION}_4V

cat ecmwf_idc_ops_header_template | sed "s,xxx,${ECUID},"  | sed "s,ggg,${ECGID},"  | sed "s,v.v,${VERSION}.${SUBVERSION}," >ecmwf_idc_ops_header
cat ecmwf_idc_ops_header ecmwf_idc_ops_body >ecmwf_idc_ops_ecgate

cat ecmwf_idc_ops_multi_header_template | sed "s,xxx,${ECUID},"  | sed "s,ggg,${ECGID},"  | sed "s,v.v,${VERSION}.${SUBVERSION}," >ecmwf_idc_ops_multi_header
cat ecmwf_idc_ops_multi_header ecmwf_idc_ops_body ecmwf_idc_ops_multi_footer >ecmwf_idc_ops_multi_ecgate


cat CONTROL_ERA_TEMPLATE | sed "s,xxx.xxx.xxx.xxx,${GATEWAY}," | sed "s,xxx@xxx,${DESTINATION}," | sed "s,PREFIX EN,PREFIX EG,"  | sed "s,v.v,${VERSION}.${SUBVERSION}," > CONTROL_ERA__GLOBALGAUSS
cat CONTROL_ERA__GLOBALGAUSS | sed "s,GAUSS 1,GAUSS 0," | sed "s,M_ETA 0,M_ETA 1," | sed "s,PREFIX EG,PREFIX EE," | sed "s,v.v,${VERSION}.${SUBVERSION},"  > CONTROL_ERA__GLOBALETA
cat CONTROL_ERA__GLOBALETA | sed "s,M_GRID 1000,M_GRID 200," | sed "s,M_RESOL 159,M_RESOL 799,"  | sed "s,M_LEFT -179000,M_LEFT -10000,"  | sed "s,M_RIGHT 180000,M_RIGHT 30000,"  | sed "s,M_LOWER -90000,M_LOWER 30000,"  | sed "s,M_UPPER 90000,M_UPPER 60000,"   | sed "s,DTIME 3,DTIME 3," | sed "s,PREFIX EE,PREFIX EH," | sed "s,v.v,${VERSION}.${SUBVERSION},"  > CONTROL_ERA__HIRES
cat CONTROL_ERA__GLOBALETA | sed "s,M_GRID 1000,M_GRID 200," | sed "s,M_RESOL 159,M_RESOL 799,"  | sed "s,M_LEFT -179000,M_LEFT 113000,"  | sed "s,M_RIGHT 180000,M_RIGHT 190000,"  | sed "s,M_LOWER -90000,M_LOWER 00000,"  | sed "s,M_UPPER 90000,M_UPPER 30000,"   | sed "s,DTIME 3,DTIME 3," | sed "s,PREFIX EE,PREFIX EH,"  | sed "s,v.v,${VERSION}.${SUBVERSION},"  > CONTROL_ERA__HAIYAN
# ERA-Interim Template
cat CONTROL_ERA__GLOBALGAUSS | sed "s,CLASS OD,CLASS EI," | sed "s,PREFIX EG,PREFIX EI," | sed "s,M_LEVEL 137,M_LEVEL 60," | sed "s,M_LEVELIST 1\/TO\/137,M_LEVELIST 1\/TO\/60,"  | sed "s,201311,201211,"  | sed "s,v.v,${VERSION}.${SUBVERSION},"  > CONTROL_ERA__EI
# ERA-Interim Template
cat CONTROL_ERA__GLOBALGAUSS | sed "s,M_TYPE AN FC FC FC FC FC AN FC FC FC FC FC AN FC FC FC FC FC AN FC FC FC FC FC,M_TYPE CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV CV," | sed "s,M_TIME 00 00 00 00 00 00 06 00 00 00 00 00 12 12 12 12 12 12 18 12 12 12 12 12,M_TIME 00 00 00 00 00 00 00 00 00 00 00 00 12 12 12 12 12 12 12 12 12 12 12 12," | sed "s,M_STEP 00 01 02 03 04 05 00 07 08 09 10 11 00 01 02 03 04 05 00 07 08 09 10 11,M_STEP 00 01 02 03 04 05 06 07 08 09 10 11 00 01 02 03 04 05 06 07 08 09 10 11,"  | sed "s,STREAM OPER,STREAM ENFO," | sed "s,M_NUMBER OFF,M_NUMBER 1,"  | sed "s,M_LEVEL 137,M_LEVEL 62," | sed "s,M_LEVELIST 1\/TO\/137,M_LEVELIST 1\/TO\/62," > CONTROL_ERA__CV

SERVER=`uname -n | cut -c 1-4`
if [[ ${SERVER} != ecgb ]] ; then
# submit via gateway from local server

echo Software copied from server `uname -n | cut -c 1-4` to ecgate

ksh upload_source V${VERSION}.${SUBVERSION}

for FILE in `ls CONTROL_ERA__*` ; do

  name=`echo $FILE | sed s,CONTROL_ERA__,,`
  cat flex_ecmwf_header $FILE flex_ecmwf_body >flex_ecmwf_$name
  ecaccess-file-put flex_ecmwf_$name flex_extract_ecgate_V${VERSION}.${SUBVERSION}/flex_ecmwf_$name

done

ecaccess-file-put CONTROL_OPS_V${VERSION}.${SUBVERSION} flex_extract_ecgate_V${VERSION}.${SUBVERSION}/CONTROL_OPS_V${VERSION}.${SUBVERSION}
ecaccess-file-put ecmwf_idc_ops_ecgate flex_extract_ecgate_V${VERSION}.${SUBVERSION}/ecmwf_idc_ops_ecgate
ecaccess-file-put ecmwf_idc_ops_multi_ecgate flex_extract_ecgate_V${VERSION}.${SUBVERSION}/ecmwf_idc_ops_multi_ecgate
ecaccess-file-mkdir scratch:ms_sms_output_V${VERSION}.${SUBVERSION}

else

for FILE in `ls CONTROL_ERA__*` ; do

  name=`echo $FILE | sed s,CONTROL_ERA__,,`
  cat flex_ecmwf_header $FILE flex_ecmwf_body >flex_ecmwf_$name
done


set +e
mkdir $SCRATCH/ms_sms_output_V${VERSION}.${SUBVERSION}
set -e

fi

echo ' '
echo !! NOTE !! NOTE !!
echo ' '
echo Scripts are now generated and uploaded but not yet submitted. 
echo Run submit_examples.ksh to submit them to ecgate.
echo ' '
echo !! NOTE !! NOTE !!
