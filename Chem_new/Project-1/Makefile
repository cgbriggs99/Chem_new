include ../delete.mk
include ../libs.mk

.PHONY:all
all: geometry_default.o

CXXFLAGS=-I $(LAPACKDIR)/include

.PHONY: clean
clean:
	$(RM) *.o

CXXFLAGS=-I /usr/local/psi4/include
LDFLAGS=-L /usr/local/psi4/lib -lint

geometry_default.o: geometry_default.cpp geometry_default.hpp geometry.hpp

geometry.hpp: ../Molecule/molecule.hpp ../Base/matrix.hpp