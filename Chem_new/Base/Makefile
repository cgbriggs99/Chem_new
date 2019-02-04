include ../libs.mk
include ../delete.mk

LDFLAGS=-L $(LAPACKDIR)/lib
LDLIBS=-lm -llapack -lgfortran -llapacke -lblas -lcblas
CXXFLAGS=-I $(LAPACKDIR)/include


.PHONY:all
all: math.o eigenvalues_default_double.o

.PHONY: clean
clean:
	$(RM) *.o

math.o: ./math.hpp

eigenvalues_default_double.o: eigenvalues_default_double.cpp eigenvalues_default.hpp eigenvalues.hpp
