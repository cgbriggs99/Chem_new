
include ../suffixes.mk

CXXFLAGS=-g3 -ggdb
EXECS=test_test$(EXEC) atom_test$(EXEC) matrix_test$(EXEC)

.PHONY: all
all: test_test atom_test matrix_test

.PHONY: run_tests
run_tests: $(EXECS)
	./test_test > ./logs/test_test.log
	./atom_test > ./logs/atom_test.log
	./matrix_test > ./logs/matrix_test.log
	
.PHONY: clean
clean:
	rm -f $(EXECS)
	rm -rf ./logs/*.log

test_test: test_test.cpp test.hpp

atom_test: atom_test.cpp ../Molecule/molecule.hpp test.hpp

matrix_test: matrix_test.cpp ../Base/matrix.hpp ../Base/matrix.cpp