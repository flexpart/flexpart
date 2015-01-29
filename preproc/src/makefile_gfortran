###############################################################################
#
# Top level Makefile for ECMWFDATA software
#
# Last modified:  August 20, 2013
#
###############################################################################

.SUFFIXES: .o .c .f .f90 .f95 .F 


.f90.o:	
	$(F90C) -c $(F90FLAGS)  $<


TARGET_ARCH =
OPT	= -g -O
DEBUG	= -g -O
LIB     = 
FC=gfortran  -m64  -fdefault-real-8 -fcray-pointer -fno-second-underscore  -ffixed-line-length-132 -fopenmp -fconvert=big-endian  -fmax-stack-var-size=100000000 
F90C=gfortran   -m64 -fdefault-real-8 -fcray-pointer -fno-second-underscore  -ffixed-line-length-132 -fopenmp  -fconvert=big-endian  -fmax-stack-var-size=100000000 
GRIB_API_INCLUDE = $(GRIB_API_FFLAGS)
GRIB_API_LIB = $(GRIB_API_LDFLAGS)
FFLAGS =  $(OPT) $(DEBUG) -I. $(GRIB_API_INCLUDE)
F90FLAGS =  $(OPT) $(DEBUG)  -I. $(GRIB_API_INCLUDE)
LDFLAGS =  $(OPT) $(DEBUG)
EMOSLIB=$(EMOS_LDFLAGS) 
LIB =  $(EMOSLIB)  $(GRIB_API_LIB) -lgcc -lm -fdefault-integer-8  -lgfortran -ljasper 

F90=$(F90C)

BINDIR  =      .

EXE	=	CONVERT2

all:	${EXE}

clean:
	rm *.o *.mod ${EXE} 

phgrreal.o: phgrreal.f90
	$(F90) -c -g -O3 -fopenmp phgrreal.f90

grphreal.o: grphreal.f90
	$(F90) -c -g -O3 -fopenmp grphreal.f90

ftrafo.o: ftrafo.f90
	$(F90) -c -g -O3 -fopenmp ftrafo.f90

CONVERT2:	phgrreal.o grphreal.o ftrafo.o rwGRIB2.o  preconvert.o
	$(F90) $(DEBUG) $(OPT) -o CONVERT2 ftrafo.o phgrreal.o grphreal.o rwGRIB2.o preconvert.o ${LIB}

#
# End of the Makefile
#
###############################################################################
