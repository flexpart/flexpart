SHELL = /bin/bash
TARGET = local
WINDS=ecmwf
#WINDS=gfs
#WINDS=fnl

FC       = gfortran
FFLAGS   =   -O2 -m64 -mcmodel=medium -fconvert=little-endian -frecord-marker=4 -I$(INCPATH)

ifeq ($(TARGET),dmz)	
  # options for ganglia
  INCPATH  = /xnilu_wrk/flex_wrk/bin64/grib_api/include
  LIBPATH1 = /xnilu_wrk/flex_wrk/bin64/grib_api/lib
  LIBPATH2 = /usr/lib/x86_64-linux-gnu/
  MAIN = FLEXPART_dmz_
endif
ifeq ($(TARGET),local)
  # local options 	
  #libs_dir=/.../flexpart/libs/
  libs_dir=/Users/ignacio/flexpart/libs/
  INCPATH  = $(libs_dir)/grib_api-1.9.9_dir/include
  LIBPATH1 = $(libs_dir)/grib_api-1.9.9_dir/lib
  LIBPATH2 = $(libs_dir)/jasper_dir/lib
  MAIN = FLEXPART_local_
endif

LDFLAGS  = $(FFLAGS) -L$(LIBPATH2) -L$(LIBPATH1) -lgrib_api_f90 -lgrib_api -lm -ljasper

MODOBJS = \
par_mod.o    com_mod.o \
conv_mod.o              hanna_mod.o \
interpol_mod.o          cmapf_mod.o \
unc_mod.o               oh_mod.o \
xmass_mod.o             flux_mod.o \
point_mod.o             outg_mod.o

OBJECTS = \
writeheader.o  writeheader_txt.o   writeheader_surf.o       assignland.o\
               part0.o \
caldate.o               partdep.o \
coordtrafo.o            psih.o \
raerod.o \
drydepokernel.o         random.o \
erf.o                   readavailable.o \
ew.o                    readcommand.o \
advance.o               readdepo.o \
releaseparticles.o      psim.o     \
FLEXPART.o              readlanduse.o \
getfields.o             init_domainfill.o\
interpol_wind.o         readoutgrid.o \
interpol_all.o          readpaths.o \
getrb.o                 readreceptors.o \
getrc.o                 readreleases.o \
getvdep.o               readspecies.o \
interpol_misslev.o      \
conccalc.o              \
concoutput.o  concoutput_surf.o          scalev.o \
pbl_profile.o           readOHfield.o\
juldate.o               timemanager.o \
interpol_vdep.o         interpol_rain.o \
partoutput.o \
hanna.o                 wetdepokernel.o \
mean.o                  wetdepo.o \
hanna_short.o           windalign.o \
hanna1.o                initialize.o \
                           calcpar_nests.o \
verttransform_nests.o   interpol_all_nests.o \
interpol_wind_nests.o   interpol_misslev_nests.o \
interpol_vdep_nests.o   interpol_rain_nests.o \
getvdep_nests.o   gridcheck_nests.o \
readwind_nests.o  \
readageclasses.o        readpartpositions.o \
calcfluxes.o            fluxoutput.o \
qvsat.o                 skplin.o \
convect43c.o               redist.o \
sort2.o                 distance.o \
centerofmass.o          plumetraj.o \
openouttraj.o           calcpv.o \
calcpv_nests.o          distance2.o \
clustering.o            interpol_wind_short.o \
interpol_wind_short_nests.o shift_field_0.o \
shift_field.o           outgrid_init.o \
openreceptors.o         boundcond_domainfill.o\
partoutput_short.o      readoutgrid_nest.o \
outgrid_init_nest.o     writeheader_nest.o writeheader_nest_surf.o \
concoutput_nest.o concoutput_surf_nest.o      wetdepokernel_nest.o \
drydepokernel_nest.o    zenithangle.o \
ohreaction.o            getvdep_nests.o \
initial_cond_calc.o     initial_cond_output.o \
dynamic_viscosity.o     get_settling.o


ifeq ($(WINDS),ecmwf)
  OBJECTS_WINDS = \
  calcpar.o          readwind.o \
  richardson.o       verttransform.o \
  obukhov.o          gridcheck.o  \
  convmix.o          calcmatrix.o 
endif

ifeq ($(WINDS),gfs)
  OBJECTS_WINDS = \
  calcpar_gfs.o          readwind_gfs.o \
  richardson_gfs.o       verttransform_gfs.o \
  obukhov_gfs.o          gridcheck_gfs.o  \
  convmix_gfs.o          calcmatrix_gfs.o
endif

$(MAIN): $(MODOBJS) $(OBJECTS)  $(OBJECTS_WINDS)
	$(FC) *.o -o $(MAIN)_$(WINDS) $(LDFLAGS)

$(OBJECTS): $(MODOBJS)

%.o: %.f90
	$(FC) -c $(FFLAGS) $<

clean:
	rm *.o *.mod

