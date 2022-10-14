.SUFFIXES: .o .f .f90
F90 = ifort
FLAGS = -g -traceback -check all -fpp

ifndef ROOTDIR 
ROOTDIR=./
endif

REALMAKEFILE=Makefile.2

MKLIBS = -L${MKLROOT}/lib -mkl 
LIBS = $(ROOTDIR)irrep_bcs.a

MODULES = comms.o init.o wave_data.o 
VECT    = main.o file_util.o 

ALL = $(MODULES) $(WAVEMD) $(VECT)

default:irvsp4

irvsp4:$(ALL)
	$(F90) $(FLAGS) -o ir2tb $(ALL) $(LIBS) $(MKLIBS)
	rm -r *.mod *.o

lib:
	(cd $(ROOTDIR)src_lib && $(MAKE) -f $(REALMAKEFILE) libs)
	cp $(ROOTDIR)src_lib/libirvsp.a .

$(MODULES):%.o:%.f90
	$(F90) -c $(FLAGS) $*.f90
$(WAVEMD):%.o:%.f90
	ifort -g -traceback -assume byterecl -c $*.f90 
$(VECT):%.o:%.f90
	$(F90) -c $(FLAGS) $*.f90

clean:
	rm -r *.mod *.o
