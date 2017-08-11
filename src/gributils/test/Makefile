
FC = gfortran

SRCDIR = ..

BINARY = testdrive

GRIBAPI = /opt/grib-api
MODOBJS = class_gribfile_mod.o
MODULES = class_gribfile_mod.mod

VPATH = ${SRCDIR}

FFLAGS = -mcmodel=medium


INCLUDES = -I${GRIBAPI}/include
LDFLAGS = -L${GRIBAPI}/lib -lgrib_api_f90 -lgrib_api

${BINARY} : ${BINARY}.o ${MODOBJS} 
	${FC} -o ${BINARY} ${BINARY}.o ${MODOBJS} ${LDFLAGS}

${BINARY}.o : ${BINARY}.f90 ${MODULES}


class_gribfile_mod.mod : class_gribfile_mod.o Makefile
	${FC} -c ../class_gribfile_mod.f90 ${FFLAGS} ${INCLUDES}

test : ${BINARY}
	( ./${BINARY} )

#------------------------------------------------------------------
clean :
	rm -r *.o *.mod ${BINARY}

%.o: %.f90
	${FC} -c ${FFLAGS} ${INCLUDES} $<


