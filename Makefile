
F90=ifort -warn all
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

genetics.o: numtypes.o error.o constants.o statistics.o genetics.f90
	$(F90) -c $^ -o $@

evolution.o: numtypes.o error.o constants.o genetics.o evolution.f90
	$(F90) -c $^ -o $@

lib: numtypes.o error.o constants.o statistics.o nonnum.o linear.o \
	int.o min.o time.o specialfunc.o root.o poly.o fourier.o genetics.o \
	evolution.o 
#	mv $(SRCDIR)/*.o .
#	mv $(SRCDIR)/*.mod .
	ar rcs libf90.a *.o

cleanall:
	rm *.o *.mod src/*~ libf90.a
