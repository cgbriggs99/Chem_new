include ../Windows/escape.mk
include ../libs.mk
include ../suffixes.mk
include ../delete.mk
include testcommand.mk


CXXFLAGS=-g3 -ggdb -I $(LAPACKDIR)/include
LDFLAGS=-L $(LAPACKDIR)/lib 
LDLIBS=-llapacke -llapack -lcblas -lblas -lgfortran
EXECS=hello$(EXEC) test_test$(EXEC) atom_test$(EXEC) matrix_test$(EXEC) eigenvalue_test$(EXEC) geometry_test$(EXEC)

.PHONY: all
all: hello test_test atom_test matrix_test eigenvalue_test geometry_test

.PHONY: run_tests
run_tests: $(EXECS)
	$(call TESTCOMM,./logs/test_memcheck.log) $(CURDIR)/test_test$(EXEC) > $(CURDIR)/logs/test_test.log
	$(call TESTCOMM,./logs/atom_memcheck.log) $(CURDIR)/atom_test$(EXEC) > $(CURDIR)/logs/atom_test.log
	$(call TESTCOMM,./logs/matrix_memcheck.log) $(CURDIR)/matrix_test$(EXEC) > $(CURDIR)/logs/matrix_test.log
	$(call TESTCOMM,./logs/eigenvalue_memcheck.log) $(CURDIR)/eigenvalue_test$(EXEC) > $(CURDIR)/logs/eigenvalue_test.log
	$(call TESTCOMM,./logs/geometry_memcheck.log) $(CURDIR)/geometry_test$(EXEC) > $(CURDIR)/logs/geometry_test.log
	$(call TESTCOMM,./logs/hello_memcheck.log) $(CURDIR)/hello$(EXEC) > $(CURDIR)/logs/hello_test.log
	
.PHONY: clean
clean:
	$(RM) $(EXECS)
ifneq ($(wildcard $(CURDIR)/logs/*.log),)
	$(RM) $(call escape,$(wildcard $(CURDIR)/logs/*.log))
endif

hello: hello.cpp

test_test: test_test.cpp test.hpp

atom_test: atom_test.cpp ../Molecule/molecule.hpp test.hpp

matrix_test: matrix_test.cpp ../Base/matrix.hpp ../Base/matrix.cpp

eigenvalue_test: eigenvalue_test.cpp ../Base/matrix.hpp ../Base/matrix.cpp ../Base/math.hpp ../Base/math.cpp ../Base/eigenvalues.hpp ../Base/eigenvalues_default.hpp ../Base/eigenvalues_default_double.cpp

geometry_test: geometry_test.cpp ../Project-1/geometry.hpp ../Project-1/geometry_default.hpp \
	../Project-1/geometry_default.cpp test.hpp ../Molecule/molecule.hpp ../Molecule/molecule.cpp \
	../Molecule/molecule_default.hpp ../Base/matrix.hpp ../Base/matrix.cpp ../Base/eigenvalues.hpp \
	../Base/eigenvalues_default.hpp ../Base/eigenvalues_default_double.cpp ../Base/math.o