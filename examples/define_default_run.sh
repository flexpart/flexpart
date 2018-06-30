#!/bin/bash
# Ignacio Pisso 2018-03-29
# generate FLEXPART variations on the default for parameter sweep hello world tests

#define input options
options_repo=~/repos/flexpart/options
#define winds
pathnames_local=~/repos/flexpart/pathnames

#copy default
cp -r  $options_repo options
#sed -i '/IOUT=/c\ IOUT=    2,' options2/COMMAND
mkdir output
cp $pathnames_local pathnames
#sed -i "s/\boptions\b/options2/g" pathnames2
#sed -i "s/\boutput\b/output2/g" pathnames2

#variations on base

# A) forward runs: 1,2,3,5
# IOUT=  1 (conc)
echo '# IOUT=  1 (conc)'
cp -r  options options1
# no change
mkdir output1
cp pathnames pathnames1
#change 2 lines in pathnames1 
sed -i "s/\boptions\b/options1/g" pathnames1
sed -i "s/\boutput\b/output1/g" pathnames1


# IOUT=    2 (mr)
cp -r  options options2
sed -i '/IOUT=/c\ IOUT=    2,' options2/COMMAND
mkdir output2
cp pathnames pathnames2
sed -i "s/\boptions\b/options2/g" pathnames2
sed -i "s/\boutput\b/output2/g" pathnames2
# IOUT=    3 (conc and mr)
cp -r  options options3
sed -i '/IOUT=/c\ IOUT=    3,' options3/COMMAND
mkdir output3
cp pathnames pathnames3
sed -i "s/\boptions\b/options3/g" pathnames3
sed -i "s/\boutput\b/output3/g" pathnames3
#../../src/FLEXPART pathnames3
#IOUT=   4 plume trajectories
# IOUT=    5 (conc +  plume trajectories)
# plume trajectories
# input
cp -r  options options5
sed -i '/IOUT=/c\ IOUT=    5,' options5/COMMAND
# output
mkdir output5
# paths
cp pathnames pathnames5
sed -i "s/\boptions\b/options5/g" pathnames5
sed -i "s/\boutput\b/output5/g" pathnames5


# launch
#nohup ../../src/FLEXPART pathnames5 &
#nohup ../../src/FLEXPART pathnames2 > nohup.out2 2>&1 &
#nohup ../../src/FLEXPART pathnames3 > nohup.out3 2>&1 &
#nohup ../../src/FLEXPART pathnames5 > nohup.out5 2>&1 &

#B) NetCDF: 9, 10, 11, 12
# IOUT=  9 (conc netCDF)
cp -r  options options9
sed -i '/IOUT=/c\ IOUT=    9,' options9/COMMAND
mkdir output9
cp pathnames pathnames9
sed -i "s/\boptions\b/options9/g" pathnames9
sed -i "s/\boutput\b/output9/g" pathnames9
# IOUT=  10 (VMR netCDF)
cp -r  options options10
sed -i '/IOUT=/c\ IOUT=    10,' options10/COMMAND
mkdir output10
cp pathnames pathnames10
sed -i "s/\boptions\b/options10/g" pathnames10
sed -i "s/\boutput\b/output10/g" pathnames10
# 11 both in netcdf
cp -r  options options11
sed -i '/IOUT=/c\ IOUT=    11,' options11/COMMAND
mkdir output11
cp pathnames pathnames11
sed -i "s/\boptions\b/options11/g" pathnames11
sed -i "s/\boutput\b/output11/g" pathnames11
# 12 
cp -r  options options12
sed -i '/IOUT=/c\ IOUT=    12,' options12/COMMAND
mkdir output12
cp pathnames pathnames12
sed -i "s/\boptions\b/options12/g" pathnames12
sed -i "s/\boutput\b/output12/g" pathnames12

#nohup ../../src/FLEXPART pathnames9 > nohup.out9 2>&1 &
#ps -e | grep FLEXPART
#FLEXPART=FLEXPART_8d70e43
#nohup $FLEXPART pathnames2 > nohup.out2 2>&1 &
#nohup $FLEXPART pathnames3 > nohup.out3 2>&1 &
#nohup $FLEXPART pathnames5 > nohup.out5 2>&1 &
#nohup $FLEXPART pathnames9 > nohup.out9 2>&1 &
#nohup $FLEXPART pathnames10 > nohup.out10 2>&1 &


#C: backward runs: bwd, bwd5, bwd_nc
#rm -r  *_bwd
#C1: bwd -- default (conc)
cp -r  options options_bwd
sed -i '/LDIRECT=/c\ LDIRECT=    -1,' options_bwd/COMMAND
sed -i '/IOUTPUTFOREACHRELEASE=/c\ IOUTPUTFOREACHRELEASE=    1,' options_bwd/COMMAND
sed -i '/IOUT=/c\ IOUT=    1,' options_bwd/COMMAND #should not be needed
mkdir output_bwd
cp pathnames pathnames_bwd
sed -i "s/\boptions\b/options_bwd/g" pathnames_bwd
sed -i "s/\boutput\b/output_bwd/g" pathnames_bwd
#C2: bwd5 -- bwd cluster 
options_new=options_bwd5
cp -r  options $options_new
sed -i '/LDIRECT=/c\ LDIRECT=    -1,' $options_new/COMMAND
sed -i '/IOUTPUTFOREACHRELEASE=/c\ IOUTPUTFOREACHRELEASE=    1,' $options_new/COMMAND
sed -i '/IOUT=/c\ IOUT=    5,' $options_new/COMMAND
output_new=output_bwd5
mkdir $output_new
pathnames_new=pathnames_bwd5
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new
#C3: bwd_nc -- bwd netCDF
suffix=_bwd_nc
options_new=options$suffix
cp -r  options $options_new
sed -i '/LDIRECT=/c\ LDIRECT=    -1,' $options_new/COMMAND
sed -i '/IOUTPUTFOREACHRELEASE=/c\ IOUTPUTFOREACHRELEASE=    1,' $options_new/COMMAND
sed -i '/IOUT=/c\ IOUT=    9,' $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new

#D: particle output
#nohup $FLEXPART pathnames_bwd > nohup.out_bwd 2>&1 &
#D1: part1: Trajectories from default (conc)
options_new=options_part1
cp -r  options $options_new
sed -i '/IPOUT=/c\ IPOUT=    1,' $options_new/COMMAND
output_new=output_part1
mkdir $output_new
pathnames_new=pathnames_part1
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new
#D2: part2 -- default + partposit in the end
options_new=options_part2
cp -r  options $options_new
sed -i '/IPOUT=/c\ IPOUT=    2,' $options_new/COMMAND
output_new=output_part2
mkdir $output_new
pathnames_new=pathnames_part2
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new
#D3: part_bwd1 --  backward trajectories 
options_new=options_part_bwd1
cp -r  options $options_new
sed -i '/LDIRECT=/c\ LDIRECT=    -1,' $options_new/COMMAND
sed -i '/IOUTPUTFOREACHRELEASE=/c\ IOUTPUTFOREACHRELEASE=    1,' $options_new/COMMAND
sed -i '/IPOUT=/c\ IPOUT=    1,' $options_new/COMMAND
output_new=output_part_bwd1
mkdir $output_new
pathnames_new=pathnames_part_bwd1
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new


#D4: part_QUASILAG: Trajectories MQUASILAG
suffix=part_QUASILAG
options_new=options$suffix
cp -r  options $options_new
sed -i '/MQUASILAG=/c\ MQUASILAG=    1,' $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new



#E: unit indices
## forward 
#E1: ind 1 2 
suffix=_ind_1_2
options_new=options$suffix
cp -r  options $options_new
sed -i '/IND_RECEPTOR=/c\ IND_RECEPTOR=    2,' $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new
#E2: ind 2 1  
suffix=_ind_2_1
#sed_in='"/IND_SOURCE=/c\ IND_SOURCE=    2,"'
options_new=options$suffix
cp -r  options $options_new
#sed -i $sed_in $options_new/COMMAND
sed -i "/IND_SOURCE=/c\ IND_SOURCE=    2," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new
#E3: ind 2 2
suffix=_ind_2_2
sed_in='"/IND_SOURCE=/c\ IND_SOURCE=    2,"'
sed_in2='"/IND_RECEPTOR=/c\ IND_RECEPTOR=    2,"'
options_new=options$suffix
cp -r  options $options_new
#sed -i $sed_in $options_new/COMMAND
#sed -i $sed_in2 $options_new/COMMAND
sed -i "/IND_SOURCE=/c\ IND_SOURCE=    2," $options_new/COMMAND
sed -i "/IND_RECEPTOR=/c\ IND_RECEPTOR=    2," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new


## backward
#E4: bwd 1 2 --   
suffix=_bwd_ind_1_2
options_new=options$suffix
cp -r  options $options_new
sed -i "/IND_RECEPTOR=/c\ IND_RECEPTOR=    2," $options_new/COMMAND
sed -i "/LDIRECT=/c\ LDIRECT=    -1," $options_new/COMMAND
sed -i "/IOUTPUTFOREACHRELEASE=/c\ IOUTPUTFOREACHRELEASE=    1," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new

#E5: bwd 2 1 -- 
suffix=_bwd_ind_2_1
options_new=options$suffix
cp -r  options $options_new
sed -i "/IND_SOURCE=/c\ IND_SOURCE=    2," $options_new/COMMAND
sed -i "/LDIRECT=/c\ LDIRECT=    -1," $options_new/COMMAND
sed -i "/IOUTPUTFOREACHRELEASE=/c\ IOUTPUTFOREACHRELEASE=    1," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new

#E6: bwd 2 2
suffix=_bwd_ind_2_2
options_new=options$suffix
cp -r  options $options_new
sed -i "/IND_SOURCE=/c\ IND_SOURCE=    2," $options_new/COMMAND
sed -i "/IND_RECEPTOR=/c\ IND_RECEPTOR=    2," $options_new/COMMAND
sed -i "/LDIRECT=/c\ LDIRECT=    -1," $options_new/COMMAND
sed -i "/IOUTPUTFOREACHRELEASE=/c\ IOUTPUTFOREACHRELEASE=    1," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new



#F different species

#F1: 
suffix=_specNO
options_new=options$suffix
cp -r  options $options_new
sed -i "/SPECNUM_REL=/c\ SPECNUM_REL=   3," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new

#F2 
suffix=_specAERO-TRACE
options_new=options$suffix
cp -r  options $options_new
sed -i "/SPECNUM_REL=/c\ SPECNUM_REL=   25," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new

#F3: 
suffix=_specCO
options_new=options$suffix
cp -r  options $options_new
sed -i "/SPECNUM_REL=/c\ SPECNUM_REL=   22," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new

#F4: 
suffix=_specBC
options_new=options$suffix
cp -r  options $options_new
sed -i "/SPECNUM_REL=/c\ SPECNUM_REL=   40," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new



#F5: 
suffix=_bwd_specNO
options_new=options$suffix
cp -r  options $options_new
sed -i "/SPECNUM_REL=/c\ SPECNUM_REL=   3," $options_new/COMMAND
sed -i "/LDIRECT=/c\ LDIRECT=    -1," $options_new/COMMAND
sed -i "/IOUTPUTFOREACHRELEASE=/c\ IOUTPUTFOREACHRELEASE=    1," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new


#G nested output
#G1
suffix=_nested
options_new=options$suffix
cp -r  options $options_new
sed -i "/NESTED_OUTPUT=/c\ NESTED_OUTPUT=   1," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new
#G2
suffix=_nested_bwd
options_new=options$suffix
cp -r  options $options_new
sed -i "/NESTED_OUTPUT=/c\ NESTED_OUTPUT=   1," $options_new/COMMAND
sed -i "/LDIRECT=/c\ LDIRECT=    -1," $options_new/COMMAND
sed -i "/IOUTPUTFOREACHRELEASE=/c\ IOUTPUTFOREACHRELEASE=    1," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new


# H  MDOMAINFILL
# H1: DOMAINFILL
suffix=_DOMAINFILL
options_new=options$suffix
cp -r  options $options_new
sed -i "/MDOMAINFILL=/c\ MDOMAINFILL=   1," $options_new/COMMAND
#sed -i "/LDIRECT=/c\ LDIRECT=    -1," $options_new/COMMAND
#sed -i "/IOUTPUTFOREACHRELEASE=/c\ IOUTPUTFOREACHRELEASE=    1," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new


# I  LINIT_COND + SURF_ONLY
# I1: Gmode
suffix=_Gmode
options_new=options$suffix
cp -r  options $options_new
sed -i "/LINIT_COND=/c\ LINIT_COND=   1," $options_new/COMMAND
sed -i "/SURF_ONLY=/c\ SURF_ONLY=   1," $options_new/COMMAND
#sed -i "/LDIRECT=/c\ LDIRECT=    -1," $options_new/COMMAND
#sed -i "/IOUTPUTFOREACHRELEASE=/c\ IOUTPUTFOREACHRELEASE=    1," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new

# I2: Gmode_bwd
suffix=_Gmode_bwd
options_new=options$suffix
cp -r  options $options_new
sed -i "/LINIT_COND=/c\ LINIT_COND=   1," $options_new/COMMAND
sed -i "/SURF_ONLY=/c\ SURF_ONLY=   1," $options_new/COMMAND
sed -i "/LDIRECT=/c\ LDIRECT=    -1," $options_new/COMMAND
sed -i "/IOUTPUTFOREACHRELEASE=/c\ IOUTPUTFOREACHRELEASE=    1," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new



# J CBLFLAG
# J 1 -- CBLFLAG fwd
suffix=_CBLFLAG
options_new=options$suffix
cp -r  options $options_new
sed -i "/CBLFLAG=/c\ CBLFLAG=   1," $options_new/COMMAND
#sed -i "/LDIRECT=/c\ LDIRECT=    -1," $options_new/COMMAND
#sed -i "/IOUTPUTFOREACHRELEASE=/c\ IOUTPUTFOREACHRELEASE=    1," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new
# J 2 -- CBLFLAG bwd
suffix=_CBLFLAG_bwd
options_new=options$suffix
cp -r  options $options_new
sed -i "/CBLFLAG=/c\ CBLFLAG=   1," $options_new/COMMAND
sed -i "/LDIRECT=/c\ LDIRECT=    -1," $options_new/COMMAND
sed -i "/IOUTPUTFOREACHRELEASE=/c\ IOUTPUTFOREACHRELEASE=    1," $options_new/COMMAND
output_new=output$suffix
mkdir $output_new
pathnames_new=pathnames$suffix
cp pathnames $pathnames_new
sed -i "s/\boptions\b/$options_new/g" $pathnames_new
sed -i "s/\boutput\b/$output_new/g" $pathnames_new





# H nested input

# need a different pathnames file

# I NCEP winds


#options_new=???  #options_part2
#cp -r  options $options_new
#???#sed -i '/IPOUT=/c\ IPOUT=    2,' $options_new/COMMAND
#output_new=???#output_part2
#mkdir $output_new
#pathnames_new=??? pathnames_part2
#cp pathnames $pathnames_new
#sed -i "s/\boptions\b/$options_new/g" $pathnames_new
#sed -i "s/\boutput\b/$output_new/g" $pathnames_new

