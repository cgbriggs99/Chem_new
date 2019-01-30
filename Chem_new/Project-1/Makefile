include ../delete.mk
include ../libs.mk

.PHONY:all
all: geometry_default.o

CXXFLAGS=-I $(LAPACKDIR)/include

.PHONY: clean
clean:
	$(RM) *.o

geometry_default.o: geometry_default.cpp geometry_default.hpp geometry.hpp

geometry.hpp: ../Molecule/molecule.hpp ../Base/matrix.hpp