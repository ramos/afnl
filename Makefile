
F90=ifort -warn all -std03
#F90=gfortran-4.3 -std=f2003 -Wall -pedantic
F90OPT=

VPATH=src
SRCDIR=$(VPATH)


all: lib

numtypes.o: numtypes.f90
	$(F90) -c $^ -o $@

constants.o: numtypes.o constants.f90
	$(F90) -c $^ -o $@

error.o: numtypes.o constants.o error.f90
	$(F90) -c $^ -o $@

nonnum.o: numtypes.o error.o constants.o nonnum.f90
	$(F90) -c $^ -o $@

linear.o: numtypes.o error.o constants.o linear.f90
	$(F90) -c $^ -o $@

int.o: numtypes.o error.o constants.o int.f90
	$(F90) -c $^ -o $@

time.o: numtypes.o error.o constants.o time.f90
	$(F90) -c $^ -o $@

specialfunc.o: numtypes.o error.o constants.o specialfunc.f90
	$(F90) -c $^ -o $@

root.o: numtypes.o error.o constants.o root.f90
	$(F90) -c $^ -o $@

min.o: numtypes.o error.o constants.o min.f90
	$(F90) -c $^ -o $@

statistics.o: numtypes.o error.o constants.o linear.o nonnum.o statistics.f90
	$(F90) -c $^ -o $@

poly.o: numtypes.o error.o constants.o poly.f90
	$(F90) -c $^ -o $@

fourier.o: numtypes.o error.o constants.o fourier.f90
	$(F90) -c $^ -o $@

minuitAPI.o: numtypes.o minuitAPI.f90
	$(F90) -c $^ -o $@

lapackAPI.o: numtypes.o lapackAPI.f90
	$(F90) -c $^ -o $@

lib: numtypes.o error.o constants.o specialfunc.o statistics.o \
	nonnum.o linear.o int.o min.o time.o root.o poly.o fourier.o \
	minuitAPI.o lapackAPI.o
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
