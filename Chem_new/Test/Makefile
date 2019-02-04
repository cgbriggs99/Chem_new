include ../Windows/escape.mk
include ../libs.mk
include ../suffixes.mk
include ../delete.mk
include testcommand.mk


CXXFLAGS=-g3 -ggdb -I $(LAPACKDIR)/include
LDFLAGS=-L $(LAPACKDIR)/lib 
LDLIBS=-llapacke -llapack -lcblas -lblas -lgfortran
EXECS=hello$(EXEC) test_test$(EXEC) atom_test$(EXEC) matrix_test$(EXEC) eigenvalue_test$(EXEC) geometry_test$(EXEC) hessian_test$(EXEC)

.PHONY: all
all: test_test atom_test matrix_test geometry_test eigenvalue_test hessian_test


.PHONY: run_tests
run_tests: $(EXECS)
	$(call TESTCOMM,./logs/test_memcheck.log) $(CURDIR)/test_test$(EXEC) > $(CURDIR)/logs/test_test.log
	$(call TESTCOMM,./logs/atom_memcheck.log) $(CURDIR)/atom_test$(EXEC) > $(CURDIR)/logs/atom_test.log
	$(call TESTCOMM,./logs/matrix_memcheck.log) $(CURDIR)/matrix_test$(EXEC) > $(CURDIR)/logs/matrix_test.log
	$(call TESTCOMM,./logs/eigenvalue_memcheck.log) $(CURDIR)/eigenvalue_test$(EXEC) > $(CURDIR)/logs/eigenvalue_test.log
	$(call TESTCOMM,./logs/geometry_memcheck.log) $(CURDIR)/geometry_test$(EXEC) > $(CURDIR)/logs/geometry_test.log
	$(call TESTCOMM,./logs/hessian_memcheck.log) $(CURDIR)/hessian_test$(EXEC) > $(CURDIR)/logs/hessian_test.log
	
.PHONY: clean
clean:
	$(RM) $(EXECS)
ifneq ($(wildcard $(CURDIR)/logs/*.log),)
	$(RM) $(call escape,$(wildcard $(CURDIR)/logs/*.log))
endif


test_test: test_test.cpp test.hpp

atom_test: atom_test.cpp ../Molecule/molecule.hpp test.hpp

matrix_test: matrix_test.cpp ../Base/matrix.hpp ../Base/matrix.cpp

eigenvalue_test: eigenvalue_test.cpp ../Base/matrix.hpp ../Base/math.hpp ../Base/eigenvalues.hpp \
	../Base/eigenvalues_default.hpp ../Base/math.o ../Base/eigenvalues_default_double.o ../Base/matrix.cpp

geometry_test: geometry_test.cpp ../Project-1/geometry_default.hpp \
	test.hpp ../Molecule/molecule.hpp ../Molecule/molecule.o ../Project-1/geometry_default.o \
	../Molecule/molecule_default.hpp ../Base/matrix.hpp ../Base/eigenvalues.hpp ../Base/math.o \
	../Base/eigenvalues_default.hpp ../Base/eigenvalues_default_double.o ../Base/matrix.cpp
	$(info Building geometry_test)
	$(CXX) $(CXXFLAGS) -o geometry_test $^ $(LDFLAGS) $(LDLIBS)
	
hessian_test: hessian_test.cpp ../Project-2/hessian_default.hpp test.hpp ../Molecule/molecule_default.hpp \
	../Base/matrix_default.hpp ../Base/eigenvalues_default.hpp ../Base/matrix.cpp ../Base/eigenvalues_default_double.o \
	../Base/math.o ../Molecule/molecule.o

