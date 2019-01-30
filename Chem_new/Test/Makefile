
include ../suffixes.mk

CXXFLAGS=-g3 -ggdb -I /usr/local/psi4/include
LDLIBS=-llapacke -L /usr/local/psi4/lib -lint -lderiv -lxc -lgg
EXECS=test_test$(EXEC) atom_test$(EXEC) matrix_test$(EXEC) eigenvalue_test$(EXEC) geometry_test$(EXEC)
VGFLAGS=--leak-check=full --show-leak-kinds=all -v

.PHONY: all
all: test_test atom_test matrix_test eigenvalue_test geometry_test

.PHONY: run_tests
run_tests: $(EXECS)
	valgrind $(VGFLAGS) --log-file="./logs/test_memcheck.log" ./test_test > ./logs/test_test.log
	valgrind $(VGFLAGS) --log-file="./logs/atom_memcheck.log" ./atom_test > ./logs/atom_test.log
	valgrind $(VGFLAGS) --log-file="./logs/matrix_memcheck.log" ./matrix_test > ./logs/matrix_test.log
	valgrind $(VGFLAGS) --log-file="./logs/eigenvalue_memcheck.log" ./eigenvalue_test > ./logs/eigenvalue_test.log
	valgrind $(VGFLAGS) --log-file="./logs/geometry_memcheck.log" ./geometry_test > ./logs/geometry_test.log
	
.PHONY: clean
clean:
	rm -f $(EXECS)
	rm -rf ./logs/*.log

test_test: test_test.cpp test.hpp

atom_test: atom_test.cpp ../Molecule/molecule.hpp test.hpp

matrix_test: matrix_test.cpp ../Base/matrix.hpp ../Base/matrix.cpp

eigenvalue_test: eigenvalue_test.cpp ../Base/matrix.hpp ../Base/matrix.cpp ../Base/math.hpp ../Base/math.cpp ../Base/eigenvalues.hpp ../Base/eigenvalues_default.hpp ../Base/eigenvalues_default_double.cpp

geometry_test: geometry_test.cpp ../Project-1/geometry.hpp ../Project-1/geometry_default.hpp \
	../Project-1/geometry_default.cpp test.hpp ../Molecule/molecule.hpp ../Molecule/molecule.cpp \
	../Molecule/molecule_default.hpp ../Base/matrix.hpp ../Base/matrix.cpp ../Base/eigenvalues.hpp \
	../Base/eigenvalues_default.hpp ../Base/eigenvalues_default_double.cpp ../Base/math.o