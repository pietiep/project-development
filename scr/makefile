# makefile
INTEL= -mkl -O3 -xW -traceback -check bounds

FC=ifort
FFLAGS= $(INTEL)# -openmp

LIB= ewpd.a
#GENLIB= ../prog/lib.a

ALL	= 							\
sil.o krytrafo.o lanczos.o l_diag.o input.o matrix.o dsyev-ex.o potential.o \
startwelle.o math_lib.o plot.o

ALLREF	= 							\
sil.f krytrafo.f lanczos.f l_diag.f input.f matrix.f dsyev-ex.f potential.f \
startwelle.f math_lib.f plot.f

EWPD	= ewpd.f $(LIB) $(ALLREF)

all:; 
	ar -xo $(LIB)
	make $(ALL)
	ar -rlc $(LIB) *.o
	rm *.o

new	: $(ALL);  
	ar -rlc $(LIB) *.o
	rm *.o

%	: $(EWPD) %.f;
	make all	
	$(FC) $(FFLAGS) -o ../prog/$* $*.f  $(LIB)

%.exe	: $(EWPD) Utilities/%.f;
	make all	
	$(FC) $(FFLAGS) -o ../prog/$*  Utilities/$*.f  $(LIB) $(GENLIB) 

