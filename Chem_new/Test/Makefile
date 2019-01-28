
include ../suffixes.mk

.PHONY:all
all: test_test atom_test

.PHONY:run_tests
run_tests: test_test$(EXEC) atom_test$(EXEC)
	test_test > logs/test_test.log
	atom_test > logs/atom_test.log

test_test: test_test.cpp test.hpp

atom_test: atom_test.cpp ../Molecule/molecule.hpp test.hpp

test.hpp:
	touch test.hpp