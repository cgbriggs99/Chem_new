.PHONY:all
all: molecule.o

.PHONY: clean
clean:
	rm -rf *.o

molecule.o: molecule.cpp molecule.hpp molecule_default.hpp

molecule.hpp: ../Base/matrix.hpp