SHELL = /bin/bash
################################################################################
#  DESCRIPTION
#    Makefile for FLEXPART. Standard (serial) and parallel (MPI) version
#
#    Dependencies are resolved in this makefile, so parallel make is
#    possible ("make -j")
# 
#  USAGE
#    Compile serial FLEXPART (ECMWF)
#      make [-j] ecmwf
#
#    Compile parallel FLEXPART (ECMWF)
#      make [-j] ecmwf-mpi
#     
#    Compile for debugging parallel FLEXPART (ECMWF)
#      make [-j] ecmwf-mpi-dbg
#
#    Compile serial FLEXPART (GFS)
#      make [-j] gfs
#
#    Compile parallel FLEXPART (GFS)
#      make [-j] gfs-mpi
#
################################################################################



## PROGRAMS
FLEXPART-ECMWF-MPI 	= FP_ecmwf_MPI
FLEXPART-ECMWF-MPI-DBG 	= DBG_FP_ecmwf_MPI
FLEXPART-ECMWF 		= FP_ecmwf_gfortran
FLEXPART-GFS 		= FP_gfs_gfortran
FLEXPART-GFS-MPI 	= FP_gfs_MPI


## COMPILERS
F90       = /usr/bin/gfortran
#F90	  = ${HOME}/gcc-4.9.1/bin/gfortran
#MPIF90    = ${HOME}/opt/bin/mpifort
#MPIF90    = mpif90.mpich
MPIF90    = mpif90.openmpi
#MPIF90    = mpifort

## OPTIMIZATION LEVEL
O_LEV = 2 # [0,1,2,3,g,s,fast]
O_LEV_DBG = 0 # [0,g]

## LIBRARIES
LIBS = -lgrib_api_f90 -lgrib_api -lm -ljasper -lnetcdff   # -fopenmp # -llapack -lnetcdf

## 1) System libraries at NILU
INCPATH1 = /usr/include
#INCPATH1 = /xnilu_wrk/flex_wrk/bin64/grib_api/include
LIBPATH1 = /xnilu_wrk/flex_wrk/bin64/grib_api/lib
#LIBPATH2 = /usr/lib/x86_64-linux-gnu

## 2) Home-made libraries
# INCPATH1  = ${HOME}/include
INCPATH2  = /homevip/flexpart/include/
LIBPATH2 = /homevip/flexpart/lib/
# LIBPATH1 = ${HOME}/lib

FFLAGS   = -I$(INCPATH1) -I$(INCPATH2) -O$(O_LEV) -g -m64 -mcmodel=medium -fconvert=little-endian -frecord-marker=4 -fmessage-length=0 -flto=jobserver -O$(O_LEV) -mtune=native -fuse-linker-plugin $(FUSER) # -march=native

DBGFLAGS = -I$(INCPATH1) -I$(INCPATH2) -O$(O_LEV_DBG) -g3 -ggdb3 -m64 -mcmodel=medium -fconvert=little-endian -frecord-marker=4 -fmessage-length=0 -flto=jobserver -O$(O_LEV_DBG) -fbacktrace -Warray-bounds  -Wall -fcheck=all $(FUSER)  #  -ffpe-trap=invalid,overflow,denormal,underflow,zero -fdump-core

LDFLAGS  = $(FFLAGS) -L$(LIBPATH1) -L$(LIBPATH2) $(LIBS)
LDDEBUG  = $(DBGFLAGS) -L$(LIBPATH1) -L$(LIBPATH2) $(LIBS)

MODOBJS = \
par_mod.o    		com_mod.o \
conv_mod.o              hanna_mod.o \
interpol_mod.o          cmapf_mod.o \
unc_mod.o               oh_mod.o \
xmass_mod.o             flux_mod.o \
point_mod.o             outg_mod.o \
	 		random_mod.o

MPI_MODOBJS = \
mpi_mod.o

## Serial versions (MPI version with same functionality and name '_mpi.f90' exists)
OBJECTS_SERIAL = \
	releaseparticles.o 	partoutput.o \
	conccalc.o \
	init_domainfill.o 	concoutput.o  \
	timemanager.o FLEXPART.o	\
	readpartpositions.o \
	partoutput_short.o		\
	concoutput_nest.o 	\
	boundcond_domainfill.o	\
	redist.o		\
	concoutput_surf.o	concoutput_surf_nest.o	\
	getfields.o

## For MPI version 
OBJECTS_MPI = releaseparticles_mpi.o partoutput_mpi.o \
	conccalc_mpi.o \
	init_domainfill_mpi.o concoutput_mpi.o 	\
	timemanager_mpi.o FLEXPART_MPI.o	\
	readpartpositions_mpi.o		\
	partoutput_short_mpi.o	 	\
	concoutput_nest_mpi.o 		\
	boundcond_domainfill_mpi.o	\
	redist_mpi.o		\
	concoutput_surf_mpi.o	concoutput_surf_nest_mpi.o	\
	getfields_mpi.o

### WINDFIELDS
## For ECMWF (serial) version:
OBJECTS_ECMWF = \
	calcpar.o          readwind.o \
	richardson.o       verttransform.o \
	obukhov.o          gridcheck.o  \
	convmix.o          calcmatrix.o


## For ECMWF MPI version:
OBJECTS_ECMWF_MPI = \
	gridcheck.o	   readwind_mpi.o \
	calcpar.o          \
	richardson.o       verttransform.o \
	obukhov.o          \
	convmix.o      calcmatrix.o

## For GFS (serial) version:
OBJECTS_GFS = \
	calcpar_gfs.o          readwind_gfs.o \
	richardson_gfs.o       verttransform_gfs.o \
	obukhov_gfs.o          gridcheck_gfs.o  \
	convmix_gfs.o          calcmatrix_gfs.o

OBJECTS = \
advance.o		initialize.o		\
writeheader.o		writeheader_txt.o  	\
writeheader_surf.o     	assignland.o\
part0.o 		gethourlyOH.o\
caldate.o               partdep.o \
coordtrafo.o            psih.o \
raerod.o 		readcommand.o 	\
drydepokernel.o         readreceptors.o \
erf.o                   readavailable.o \
ew.o			readreleases.o  \
readdepo.o \
psim.o			outgrid_init.o 	\
outgrid_init_nest.o   	\
photo_O1D.o 		readlanduse.o \
interpol_wind.o         readoutgrid.o \
interpol_all.o          readpaths.o \
getrb.o                 \
getrc.o                 \
getvdep.o               readspecies.o \
interpol_misslev.o      \
scalev.o \
pbl_profile.o           readOHfield.o\
juldate.o               \
interpol_vdep.o         interpol_rain.o \
hanna.o                 wetdepokernel.o \
mean.o                  wetdepo.o \
hanna_short.o           windalign.o \
hanna1.o                	\
                        gridcheck_nests.o \
readwind_nests.o        calcpar_nests.o \
verttransform_nests.o   interpol_all_nests.o \
interpol_wind_nests.o   interpol_misslev_nests.o \
interpol_vdep_nests.o   interpol_rain_nests.o \
readageclasses.o         \
calcfluxes.o            fluxoutput.o \
qvsat.o                 skplin.o \
convect43c.o            \
sort2.o                 distance.o \
centerofmass.o          plumetraj.o \
openouttraj.o           calcpv.o \
calcpv_nests.o          distance2.o \
clustering.o            interpol_wind_short.o \
interpol_wind_short_nests.o shift_field_0.o \
shift_field.o             \
openreceptors.o         \
readoutgrid_nest.o \
writeheader_nest.o writeheader_nest_surf.o \
wetdepokernel_nest.o \
drydepokernel_nest.o    zenithangle.o \
ohreaction.o            getvdep_nests.o \
initial_cond_calc.o     initial_cond_output.o \
dynamic_viscosity.o     get_settling.o	\
initialize_cbl_vel.o	re_initialize_particle.o \
cbl.o			netcdf_output_mod.o

%.o: %.mod

ecmwf: $(FLEXPART-ECMWF)
ecmwf: FC := $(F90)

ecmwf-mpi: $(FLEXPART-ECMWF-MPI)
ecmwf-mpi: FC := $(MPIF90)

ecmwf-mpi-dbg: $(FLEXPART-ECMWF-MPI-DBG)
ecmwf-mpi-dbg: FFLAGS := $(DBGFLAGS) 
ecmwf-mpi-dbg: LDFLAGS:= $(LDDEBUG) 
ecmwf-mpi-dbg: FC := $(MPIF90)

gfs: $(FLEXPART-GFS)
gfs: FC := $(F90)
gfs-mpi: $(FLEXPART-GFS-MPI)
gfs-mpi: FC := $(MPIF90)

#all: $(FLEXPART-ECMWF)
#all: $(FLEXPART-ECMWF-MPI)

$(FLEXPART-ECMWF): $(MODOBJS) $(OBJECTS) $(OBJECTS_SERIAL) $(OBJECTS_ECMWF)
	+$(FC) -o $@ $(MODOBJS) $(OBJECTS) $(OBJECTS_SERIAL) $(OBJECTS_ECMWF) $(LDFLAGS)

$(FLEXPART-ECMWF-MPI): $(MODOBJS) $(MPI_MODOBJS) $(OBJECTS) $(OBJECTS_MPI) $(OBJECTS_ECMWF_MPI)
	+$(FC) -o $@ $(MODOBJS) $(MPI_MODOBJS) $(OBJECTS) $(OBJECTS_MPI) \
	$(OBJECTS_ECMWF_MPI) $(LDFLAGS)
#	+$(FC) -o $@ *.o $(LDFLAGS)

$(FLEXPART-ECMWF-MPI-DBG): $(MODOBJS) $(MPI_MODOBJS) $(OBJECTS) $(OBJECTS_MPI) \
	$(OBJECTS_ECMWF_MPI)
	+$(FC) -o $@ $(MODOBJS) $(MPI_MODOBJS) $(OBJECTS) $(OBJECTS_MPI) \
	$(OBJECTS_ECMWF_MPI) $(LDFLAGS)

$(FLEXPART-GFS): $(MODOBJS) $(OBJECTS) $(OBJECTS_SERIAL) $(OBJECTS_GFS)
	+$(FC) -o $@ $(MODOBJS) $(OBJECTS) $(OBJECTS_SERIAL) $(OBJECTS_GFS) $(LDFLAGS)

$(FLEXPART-GFS-MPI): $(MODOBJS) $(MPI_MODOBJS) $(OBJECTS) $(OBJECTS_MPI) $(OBJECTS_GFS)
	+$(FC) -o $@ $(MODOBJS) $(MPI_MODOBJS) $(OBJECTS) $(OBJECTS_MPI) \
	$(OBJECTS_GFS) $(LDFLAGS)

%.o: %.f90
	+$(FC) -c $(FFLAGS) $<

clean:
	\rm -f *.o *.mod

cleanall:
	\rm -f *.o *.mod $(FLEXPART-ECMWF-MPI) $(FLEXPART-ECMWF-MPI-DBG) $(FLEXPART-ECMWF) \
	$(FLEXPART-GFS-MPI) $(FLEXPART-GFS)

.SUFFIXES = $(SUFFIXES) .f90

## DEPENDENCIES
advance.o: cmapf_mod.o com_mod.o hanna_mod.o interpol_mod.o par_mod.o \
	point_mod.o random_mod.o
assignland.o: com_mod.o par_mod.o
boundcond_domainfill.o: com_mod.o par_mod.o point_mod.o random_mod.o
boundcond_domainfill_mpi.o: com_mod.o mpi_mod.o par_mod.o point_mod.o \
	random_mod.o
calcfluxes.o: com_mod.o flux_mod.o outg_mod.o par_mod.o
calcmatrix.o: com_mod.o conv_mod.o par_mod.o
calcmatrix_gfs.o: com_mod.o conv_mod.o par_mod.o
calcpar.o: com_mod.o par_mod.o
calcpar_gfs.o: com_mod.o par_mod.o
calcpar_nests.o: com_mod.o par_mod.o
calcpv.o: com_mod.o par_mod.o
calcpv_nests.o: com_mod.o par_mod.o
caldate.o: par_mod.o
cbl.o: com_mod.o par_mod.o
centerofmass.o: par_mod.o
clustering.o: par_mod.o
cmapf_mod.o: par_mod.o
com_mod.o: par_mod.o
conccalc.o: com_mod.o outg_mod.o par_mod.o unc_mod.o
conccalc_mpi.o: com_mod.o mpi_mod.o outg_mod.o par_mod.o unc_mod.o
concoutput.o: com_mod.o outg_mod.o par_mod.o point_mod.o unc_mod.o
concoutput_mpi.o: com_mod.o mpi_mod.o outg_mod.o par_mod.o point_mod.o \
	unc_mod.o
concoutput_nest.o: com_mod.o outg_mod.o par_mod.o point_mod.o unc_mod.o
concoutput_nest_mpi.o: com_mod.o mpi_mod.o outg_mod.o par_mod.o point_mod.o \
	unc_mod.o
concoutput_surf.o: com_mod.o outg_mod.o par_mod.o point_mod.o unc_mod.o
concoutput_surf_mpi.o: com_mod.o mpi_mod.o outg_mod.o par_mod.o point_mod.o \
	unc_mod.o
concoutput_surf_nest.o: com_mod.o outg_mod.o par_mod.o point_mod.o unc_mod.o
concoutput_surf_nest_mpi.o: com_mod.o mpi_mod.o outg_mod.o par_mod.o \
	point_mod.o unc_mod.o
conv_mod.o: par_mod.o
convect43c.o: conv_mod.o par_mod.o
convmix.o: com_mod.o conv_mod.o flux_mod.o par_mod.o
convmix_gfs.o: com_mod.o conv_mod.o par_mod.o
coordtrafo.o: com_mod.o par_mod.o point_mod.o
distance.o: par_mod.o
distance2.o: par_mod.o
drydepokernel.o: com_mod.o par_mod.o unc_mod.o
drydepokernel_nest.o: com_mod.o par_mod.o unc_mod.o
erf.o: par_mod.o
FLEXPART.o: com_mod.o conv_mod.o par_mod.o point_mod.o random_mod.o netcdf_output_mod.o
FLEXPART_MPI.o: com_mod.o conv_mod.o mpi_mod.o par_mod.o point_mod.o \
	random_mod.o netcdf_output_mod.o
fluxoutput.o: com_mod.o flux_mod.o outg_mod.o par_mod.o
get_settling.o: com_mod.o par_mod.o
getfields.o: com_mod.o par_mod.o
getfields_mpi.o: com_mod.o par_mod.o mpi_mod.o
gethourlyOH.o: com_mod.o oh_mod.o par_mod.o
getrb.o: par_mod.o
getrc.o: com_mod.o par_mod.o
getvdep.o: com_mod.o par_mod.o
getvdep_nests.o: com_mod.o par_mod.o
gridcheck.o: cmapf_mod.o com_mod.o conv_mod.o par_mod.o
gridcheck_emos.o: com_mod.o conv_mod.o par_mod.o
gridcheck_fnl.o: cmapf_mod.o com_mod.o conv_mod.o par_mod.o
gridcheck_gfs.o: cmapf_mod.o com_mod.o conv_mod.o par_mod.o
gridcheck_gfs_emos.o: com_mod.o conv_mod.o par_mod.o
gridcheck_nests.o: com_mod.o par_mod.o
gridcheck_nests_emos.o: com_mod.o par_mod.o
gridcheck_orig_ecmwf.o: cmapf_mod.o com_mod.o conv_mod.o par_mod.o
hanna.o: com_mod.o hanna_mod.o par_mod.o
hanna1.o: com_mod.o hanna_mod.o par_mod.o
hanna_short.o: com_mod.o hanna_mod.o par_mod.o
init_domainfill.o: com_mod.o par_mod.o point_mod.o random_mod.o
init_domainfill_mpi.o: com_mod.o mpi_mod.o par_mod.o point_mod.o \
	random_mod.o
initial_cond_calc.o: com_mod.o outg_mod.o par_mod.o unc_mod.o
initial_cond_output.o: com_mod.o outg_mod.o par_mod.o point_mod.o unc_mod.o
initialize.o: com_mod.o hanna_mod.o interpol_mod.o par_mod.o random_mod.o
initialize_cbl_vel.o: com_mod.o par_mod.o random_mod.o
interpol_all.o: com_mod.o hanna_mod.o interpol_mod.o par_mod.o
interpol_all_nests.o: com_mod.o hanna_mod.o interpol_mod.o par_mod.o
interpol_misslev.o: com_mod.o hanna_mod.o interpol_mod.o par_mod.o
interpol_misslev_nests.o: com_mod.o hanna_mod.o interpol_mod.o par_mod.o
interpol_mod.o: par_mod.o
interpol_rain.o: par_mod.o
interpol_rain_nests.o: par_mod.o
interpol_vdep.o: com_mod.o interpol_mod.o par_mod.o
interpol_vdep_nests.o: com_mod.o interpol_mod.o par_mod.o
interpol_wind.o: com_mod.o interpol_mod.o par_mod.o
interpol_wind_nests.o: com_mod.o interpol_mod.o par_mod.o
interpol_wind_short.o: com_mod.o interpol_mod.o par_mod.o
interpol_wind_short_nests.o: com_mod.o interpol_mod.o par_mod.o
juldate.o: par_mod.o
mpi_mod.o: com_mod.o par_mod.o unc_mod.o
netcdf_output_mod.o: com_mod.o outg_mod.o par_mod.o point_mod.o unc_mod.o
obukhov.o: par_mod.o
obukhov_gfs.o: par_mod.o
ohreaction.o: com_mod.o oh_mod.o par_mod.o
openouttraj.o: com_mod.o par_mod.o point_mod.o
openreceptors.o: com_mod.o par_mod.o
outgrid_init.o: com_mod.o flux_mod.o oh_mod.o outg_mod.o par_mod.o unc_mod.o
outgrid_init_nest.o: com_mod.o outg_mod.o par_mod.o unc_mod.o
part0.o: par_mod.o
partdep.o: par_mod.o
partoutput.o: com_mod.o par_mod.o
partoutput_mpi.o: com_mod.o mpi_mod.o par_mod.o
partoutput_short.o: com_mod.o par_mod.o
partoutput_short_mpi.o: com_mod.o mpi_mod.o par_mod.o
pbl_profile.o: par_mod.o
plumetraj.o: com_mod.o par_mod.o point_mod.o
psih.o: par_mod.o
psim.o: par_mod.o
raerod.o: par_mod.o
re_initialize_particle.o: com_mod.o par_mod.o
readageclasses.o: com_mod.o par_mod.o
readavailable.o: com_mod.o par_mod.o
readcommand.o: com_mod.o par_mod.o
readdepo.o: com_mod.o par_mod.o
readlanduse.o: com_mod.o par_mod.o
readlanduse_int1.o: com_mod.o par_mod.o
readOHfield.o: com_mod.o oh_mod.o par_mod.o
readoutgrid.o: com_mod.o outg_mod.o par_mod.o
readoutgrid_nest.o: com_mod.o outg_mod.o par_mod.o
readpartpositions.o: com_mod.o par_mod.o random_mod.o
readpartpositions_mpi.o: com_mod.o mpi_mod.o par_mod.o random_mod.o
readpaths.o: com_mod.o par_mod.o
readreceptors.o: com_mod.o par_mod.o
readreleases.o: com_mod.o par_mod.o point_mod.o xmass_mod.o
readspecies.o: com_mod.o par_mod.o
readwind.o: com_mod.o par_mod.o
readwind_emos.o: com_mod.o par_mod.o
readwind_gfs.o: com_mod.o par_mod.o
readwind_gfs_emos.o: com_mod.o par_mod.o
readwind_mpi.o: com_mod.o mpi_mod.o par_mod.o
readwind_nests.o: com_mod.o par_mod.o
readwind_nests_emos.o: com_mod.o par_mod.o
redist.o: com_mod.o conv_mod.o par_mod.o random_mod.o
redist_mpi.o: com_mod.o conv_mod.o mpi_mod.o par_mod.o random_mod.o
releaseparticles.o: com_mod.o par_mod.o point_mod.o random_mod.o \
	xmass_mod.o
releaseparticles_mpi.o: com_mod.o mpi_mod.o par_mod.o point_mod.o \
	random_mod.o xmass_mod.o
richardson.o: par_mod.o
richardson_gfs.o: par_mod.o
scalev.o: par_mod.o
shift_field.o: par_mod.o
shift_field_0.o: par_mod.o
timemanager.o: com_mod.o flux_mod.o netcdf_output_mod.o oh_mod.o outg_mod.o \
	par_mod.o point_mod.o unc_mod.o xmass_mod.o
timemanager_mpi.o: com_mod.o flux_mod.o mpi_mod.o oh_mod.o outg_mod.o \
	par_mod.o point_mod.o unc_mod.o xmass_mod.o
verttransform.o: cmapf_mod.o com_mod.o par_mod.o
verttransform_gfs.o: cmapf_mod.o com_mod.o par_mod.o
verttransform_nests.o: com_mod.o par_mod.o
wetdepo.o: com_mod.o par_mod.o point_mod.o
wetdepokernel.o: com_mod.o par_mod.o unc_mod.o
wetdepokernel_nest.o: com_mod.o par_mod.o unc_mod.o
writeheader.o: com_mod.o outg_mod.o par_mod.o point_mod.o
writeheader_nest.o: com_mod.o outg_mod.o par_mod.o point_mod.o
writeheader_nest_surf.o: com_mod.o outg_mod.o par_mod.o point_mod.o
writeheader_surf.o: com_mod.o outg_mod.o par_mod.o point_mod.o
writeheader_txt.o: com_mod.o outg_mod.o par_mod.o point_mod.o
zenithangle.o: par_mod.o
