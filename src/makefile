SHELL = /bin/bash
################################################################################
#  DESCRIPTION
#    Makefile for FLEXPART. Standard (serial) and parallel (MPI) version
#
#    Dependencies are resolved in this makefile, so parallel make is
#    possible ("make -j")
#
#    At NILU we have installed gcc-4.9.1 and libraries under user /homevip/flexpart
#    ("ROOT_DIR")
#    To use gfortran version 4.9, add "gcc=4.9" to the make command, e.g.
#       'make -j ecmwf gcc=4.9',
#    also set environment variable LD_LIBRARY_PATH to point to compiler libraries
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

ifeq ($(gcc), 4.9)
# Compiled libraries under users ~flexpart, gfortran v4.9
	ROOT_DIR = /homevip/flexpart/
#	ROOT_DIR = /homevip/espen/

	F90	  = ${ROOT_DIR}/gcc-4.9.1/bin/gfortran
	MPIF90    = ${ROOT_DIR}/bin/mpifort

	INCPATH1  = ${ROOT_DIR}/gcc-4.9.1/include	
	INCPATH2  = ${ROOT_DIR}/include
	LIBPATH1 = ${ROOT_DIR}/lib

else
# Default: System libraries at NILU, gfortran v4.6
	F90       = /usr/bin/gfortran
	MPIF90    = /usr/bin/mpif90.openmpi

	INCPATH1 = /xnilu_wrk/flex_wrk/bin64/grib_api/include
	INCPATH2 = /usr/include
	LIBPATH1 = /xnilu_wrk/flex_wrk/bin64/grib_api/lib
endif


## OPTIMIZATION LEVEL
O_LEV = 2 # [0,1,2,3,g,s,fast]
O_LEV_DBG = g # [0,g]

## LIBRARIES
LIBS = -lgrib_api_f90 -lgrib_api -lm -ljasper -lnetcdff # -fopenmp

FFLAGS   = -I$(INCPATH1) -I$(INCPATH2) -O$(O_LEV) -g -m64 -mcmodel=medium -fconvert=little-endian -frecord-marker=4 -fmessage-length=0 -flto=jobserver -O$(O_LEV) $(FUSER)  #-Warray-bounds -fcheck=all # -march=native

DBGFLAGS = -I$(INCPATH1) -I$(INCPATH2) -O$(O_LEV_DBG) -g3 -ggdb3 -m64 -mcmodel=medium -fconvert=little-endian -frecord-marker=4 -fmessage-length=0 -flto=jobserver -O$(O_LEV_DBG) -fbacktrace   -Wall  -fdump-core $(FUSER)  #  -ffpe-trap=invalid,overflow,denormal,underflow,zero  -Warray-bounds -fcheck=all

LDFLAGS  = $(FFLAGS) -L$(LIBPATH1) $(LIBS) #-L$(LIBPATH2)
LDDEBUG  = $(DBGFLAGS) -L$(LIBPATH1) $(LIBS) #-L$(LIBPATH2)

MODOBJS = \
par_mod.o    		com_mod.o \
conv_mod.o              hanna_mod.o \
interpol_mod.o          cmapf_mod.o \
unc_mod.o               oh_mod.o \
xmass_mod.o             flux_mod.o \
point_mod.o             outg_mod.o \
mean_mod.o	 	random_mod.o

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
	convmix.o          calcmatrix.o \
	ecmwf_mod.o


## For ECMWF MPI version:
OBJECTS_ECMWF_MPI = \
	gridcheck.o	   readwind_mpi.o \
	calcpar.o          \
	richardson.o       verttransform.o \
	obukhov.o          \
	convmix.o	   calcmatrix.o \
	ecmwf_mod.o

## For GFS (serial) version:
OBJECTS_GFS = \
	calcpar_gfs.o          readwind_gfs.o \
	richardson_gfs.o       verttransform_gfs.o \
	obukhov_gfs.o          gridcheck_gfs.o  \
	convmix_gfs.o          calcmatrix_gfs.o \
	gfs_mod.o

OBJECTS = \
advance.o		initialize.o		\
writeheader.o		writeheader_txt.o  	\
writeprecip.o \
writeheader_surf.o     	assignland.o\
part0.o 		gethourlyOH.o\
caldate.o               partdep.o \
coordtrafo.o            psih.o \
raerod.o 		readcommand.o 	\
drydepokernel.o         readreceptors.o \
erf.o                   readavailable.o \
ew.o			readreleases.o  \
readdepo.o              get_vdep_prob.o   \
get_wetscav.o	\
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
                  wetdepo.o \
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

## This allows for switching between ECMWF/GFS without editing source code.
wind_mod = ecmwf_mod.o # default wind
# ifeq ($(MAKECMDGOALS),ecmwf)
# wind_mod = ecmwf_mod.o
# endif
# ifeq ($(MAKECMDGOALS),ecmwf-mpi)
# wind_mod = ecmwf_mod.o
# endif
# ifeq ($(MAKECMDGOALS),ecmwf-mpi-dbg)
# wind_mod = ecmwf_mod.o
# endif

ifeq ($(MAKECMDGOALS),gfs)
wind_mod = gfs_mod.o
endif
ifeq ($(MAKECMDGOALS),gfs-mpi)
wind_mod = gfs_mod.o
endif

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
get_vdep_prob.o: cmapf_mod.o com_mod.o hanna_mod.o interpol_mod.o par_mod.o \
	point_mod.o random_mod.o
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
concoutput.o: com_mod.o outg_mod.o par_mod.o point_mod.o unc_mod.o mean_mod.o
concoutput_mpi.o: com_mod.o mpi_mod.o outg_mod.o par_mod.o point_mod.o \
	unc_mod.o mean_mod.o
concoutput_nest.o: com_mod.o outg_mod.o par_mod.o point_mod.o unc_mod.o mean_mod.o
concoutput_nest_mpi.o: com_mod.o mpi_mod.o outg_mod.o par_mod.o point_mod.o \
	unc_mod.o mean_mod.o
concoutput_surf.o: com_mod.o outg_mod.o par_mod.o point_mod.o unc_mod.o mean_mod.o
concoutput_surf_mpi.o: com_mod.o mpi_mod.o outg_mod.o par_mod.o point_mod.o \
	unc_mod.o mean_mod.o
concoutput_surf_nest.o: com_mod.o outg_mod.o par_mod.o point_mod.o unc_mod.o mean_mod.o
concoutput_surf_nest_mpi.o: com_mod.o mpi_mod.o outg_mod.o par_mod.o \
	point_mod.o unc_mod.o mean_mod.o
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
mean_mod.o: par_mod.o
mpi_mod.o: com_mod.o par_mod.o unc_mod.o
netcdf_output_mod.o: com_mod.o outg_mod.o par_mod.o point_mod.o unc_mod.o mean_mod.o
obukhov.o: par_mod.o
obukhov_gfs.o: par_mod.o
ohreaction.o: com_mod.o oh_mod.o par_mod.o
openouttraj.o: com_mod.o par_mod.o point_mod.o
openreceptors.o: com_mod.o par_mod.o
outg_mod.o: par_mod.o
outgrid_init.o: com_mod.o flux_mod.o oh_mod.o outg_mod.o par_mod.o unc_mod.o
outgrid_init_nest.o: com_mod.o outg_mod.o par_mod.o unc_mod.o
par_mod.o : $(wind_mod)
part0.o: par_mod.o
partdep.o: par_mod.o
partoutput.o: com_mod.o par_mod.o
partoutput_mpi.o: com_mod.o mpi_mod.o par_mod.o
partoutput_short.o: com_mod.o par_mod.o
partoutput_short_mpi.o: com_mod.o mpi_mod.o par_mod.o
pbl_profile.o: par_mod.o
plumetraj.o: com_mod.o par_mod.o point_mod.o mean_mod.o
psih.o: par_mod.o
psim.o: par_mod.o
raerod.o: par_mod.o
re_initialize_particle.o: com_mod.o par_mod.o
readageclasses.o: com_mod.o par_mod.o
readavailable.o: com_mod.o par_mod.o
readcommand.o: com_mod.o par_mod.o
readdepo.o: com_mod.o par_mod.o
readlanduse.o: com_mod.o par_mod.o
#readlanduse_int1.o: com_mod.o par_mod.o
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
unc_mod.o: par_mod.o
verttransform.o: cmapf_mod.o com_mod.o par_mod.o 
verttransform_gfs.o: cmapf_mod.o com_mod.o par_mod.o
verttransform_nests.o: com_mod.o par_mod.o
get_wetscav.o: com_mod.o par_mod.o point_mod.o
wetdepo.o: com_mod.o par_mod.o point_mod.o
wetdepokernel.o: com_mod.o par_mod.o unc_mod.o
wetdepokernel_nest.o: com_mod.o par_mod.o unc_mod.o
writeheader.o: com_mod.o outg_mod.o par_mod.o point_mod.o
writeprecip.o: com_mod.o par_mod.o point_mod.o
writeheader_nest.o: com_mod.o outg_mod.o par_mod.o point_mod.o
writeheader_nest_surf.o: com_mod.o outg_mod.o par_mod.o point_mod.o
writeheader_surf.o: com_mod.o outg_mod.o par_mod.o point_mod.o
writeheader_txt.o: com_mod.o outg_mod.o par_mod.o point_mod.o
zenithangle.o: par_mod.o
