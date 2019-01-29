.PHONY:all
all: geometry_default.o

.PHONY: clean
clean:
	rm -rf *.o

geometry_default.o: geometry_default.cpp geometry_default.hpp geometry.hpp

geometry.hpp: ../Molecule/molecule.hpp ../Base/matrix.hpp