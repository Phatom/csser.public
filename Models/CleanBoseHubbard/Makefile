BASEPATH=../../
SSEPATH=-I $(BASEPATH)/include -I $(BASEPATH)/core -I $(BASEPATH)/modes -I $(BASEPATH)/observables -I $(BASEPATH)/simplex/nrlib

CPP=mpicxx 
FORTRAN=
FLAGS=
CFLAGS= -O3 -fopenmp -lm -std=gnu++0x -lgsl $(SSEPATH) -I ./ 
LDFLAGS = -fopenmp -lgsl -lgslcblas -L$(BASEPATH) -lsse 
F_SOURCES= 
DEBUG=
CPP_SOURCES= ./*.cpp  
OBJECTS=*.o
EXECUTABLE =  EX

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CPP) $(OBJECTS) $(LDFLAGS) -o $@

$(OBJECTS): $(CPP_SOURCES)
	$(CPP) $(CFLAGS) -c $(DEBUG) $(CPP_SOURCES)
#	$(FORTRAN) $(FFLAGS) -c $(F_SOURCES)

clean: 
	rm -f *.o
