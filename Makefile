
F90=ifort -warn all
F90OPT=

VPATH=src
SRCDIR=$(VPATH)


all: lib

numtypes.o: numtypes.f90
	$(F90) -c $^ -o $@

%.o: numtypes.o %.f90
	$(F90) -c $^ -o $@

lib: numtypes.o $(patsubst %.f90,%.o,$(wildcard $(SRCDIR)/*.f90))
	mv $(SRCDIR)/*.o .
#	mv $(SRCDIR)/*.mod .
	ar rcs libf90.a *.o

cleanall:
	rm *.o *.mod src/*~ libf90.a
