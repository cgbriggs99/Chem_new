
LDFLAGS=-lm -llapack -lgfortran -llapacke -lblas -lcblas

.PHONY:all
all: math.o eigenvalues_default_double.o

.PHONY: clean
clean:
	rm -rf *.o

math.o: ./math.hpp

eigenvalues_default_double.o: eigenvalues_default_double.cpp eigenvalues_default.hpp eigenvalues.hpp
