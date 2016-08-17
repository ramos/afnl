
#F90=ifort -warn all -std03
#F90=gfortran6 -std=f2008 -ftree-vectorize -ftree-vectorizer-verbose=2 -march=native -fopt-info -fopt-info-optimized -O3
F90=gfortran6 -std=f2008 -Wall -pedantic -march=native -funroll-loops -O0 -finline-limit=60 -fbacktrace -ffpe-trap=underflow,overflow,denormal
F90OPT=

VPATH=src
SRCDIR=$(VPATH)


all: lib

numtypes.o: numtypes.f90
	$(F90) -c $< -o $@

constants.o: constants.f90 numtypes.o
	$(F90) -c $< -o $@

error.o: error.f90 numtypes.o constants.o
	$(F90) -c $< -o $@

nonnum.o: nonnum.f90 numtypes.o error.o constants.o time.o
	$(F90) -c $< -o $@

linear.o: linear.f90 numtypes.o error.o constants.o
	$(F90) -c $< -o $@

int.o: int.f90 numtypes.o error.o constants.o
	$(F90) -c $< -o $@

time.o: time.f90 numtypes.o error.o constants.o
	$(F90) -c $< -o $@

specialfunc.o: specialfunc.f90 numtypes.o error.o constants.o
	$(F90) -c $< -o $@

root.o: root.f90 numtypes.o error.o constants.o
	$(F90) -c $< -o $@

min.o: min.f90 numtypes.o error.o constants.o
	$(F90) -c $< -o $@

statistics.o: statistics.f90 numtypes.o error.o constants.o linear.o nonnum.o
	$(F90) -c $< -o $@

poly.o: poly.f90 numtypes.o error.o constants.o
	$(F90) -c $< -o $@

fourier.o: fourier.f90 numtypes.o error.o constants.o
	$(F90) -c $< -o $@

minuitAPI.o: minuitAPI.f90 numtypes.o 
	$(F90) -c $< -o $@

lapackAPI.o: lapackAPI.f90 numtypes.o 
	$(F90) -c $< -o $@

bdio.o: bdio.f90
	$(F90) -c $< -o $@

mixmax.o: mixmax.f90
	$(F90) -c $< -o $@

ranlux.o: ranlux.f90
	$(F90) -c $< -o $@

ranlux48.o: ranlux48.f90
	$(F90) -c $< -o $@

random.o: random.f90 numtypes.o ranlux.o mixmax.o 
	$(F90) -c $< -o $@

lib: numtypes.o error.o constants.o specialfunc.o mixmax.o ranlux.o \
	ranlux48.o random.o statistics.o \
	nonnum.o linear.o int.o min.o time.o root.o poly.o fourier.o \
	minuitAPI.o lapackAPI.o bdio.o
#	mv $(SRCDIR)/*.o .
#	mv $(SRCDIR)/*.mod .
	ar rcs libafnl.a *.o

lib-no-API: numtypes.o error.o constants.o specialfunc.o statistics.o \
	nonnum.o linear.o int.o min.o time.o root.o poly.o fourier.o 
#	mv $(SRCDIR)/*.o .
#	mv $(SRCDIR)/*.mod .
	ar rcs libafnl.a *.o

cleanall:
	rm *.o *.mod src/*~ libafnl.a
