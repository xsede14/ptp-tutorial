CC = mpicc
CFLAGS = -g -fopenmp
FC = mpif90
FFLAGS = -g -fopenmp

# Flags for gprof,gcov builds
#CFLAGS = -g -pg -ftest-coverage -fprofile-arcs
#FFLAGS = -g -pg -ftest-coverage -fprofile-arcs

# determine whether mpicc invokes GNU or Intel compilers and set LIB accordingly
COMPILER=$(shell $(CC) --version | head -1 | cut -d ' ' -f 1)
$(info Note: mpicc appears to invoke $(COMPILER))
ifeq ($(COMPILER),gcc)           # Set libraries for GNU compilers
    LIB = -lgfortran
else ifeq ($(COMPILER),icc)      # Set libraries for Intel compilers
    LIB = -lifcore -limf -ldl
else
    $(error The compiler $(COMPILER) does not appear to be GCC or ICC.  Please fix the Makefile)
endif

# include file used by TAU; empty usually
include eclipse.inc

SOURCES = calc.c copy.c diag.c init.c main.c time.c tstep.f90 worker.c dump.c
OBJECTS = calc.o copy.o diag.o init.o main.o time.o tstep.o worker.o dump.o
HEADERS = decs.h

all: shallow

shallow: $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) -o shallow $(OBJECTS) -lm $(LIB)

tags: $(SOURCES) $(HEADERS)
	ctags -w $(SOURCES) $(HEADERS)

clean:
	-rm -f shallow $(OBJECTS) core

calc.o: calc.c decs.h 
copy.o: copy.c decs.h 
diag.o: diag.c decs.h 
init.o: init.c decs.h 
main.o: main.c decs.h 
time.o: time.c decs.h 

tstep.o: tstep.f90
	$(FC) $(FFLAGS)  -c -o $@ $<

worker.o: worker.c decs.h
dump.o: dump.c decs.h
