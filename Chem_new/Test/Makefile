include ../Windows/escape.mk
include ../libs.mk
include ../suffixes.mk
include ../delete.mk
include testcommand.mk


CXXFLAGS=-g3 -ggdb -I $(LAPACKDIR)/include
CFLAGS=-g3 -ggdb -I $(LAPACKDIR)/include
LDFLAGS=-L $(LAPACKDIR)/lib 
LDLIBS=-llapacke -llapack -lcblas -lblas -lgfortran
EXECS=test_test$(EXEC) atom_test$(EXEC) matrix_test$(EXEC) matrix_arithmetic_test$(EXEC) \
	eigenvalue_test$(EXEC) geometry_test$(EXEC) hessian_test$(EXEC) scf_test$(EXEC)

.PHONY: all
all: test_test atom_test matrix_test matrix_arithmetic_test geometry_test eigenvalue_test hessian_test scf_test


.PHONY: run_tests
run_tests: $(EXECS)
	$(call TESTCOMM,./logs/test_memcheck.log) $(call escape,$(CURDIR)/test_test$(EXEC) > $(CURDIR)/logs/test_test.log)
	$(call TESTCOMM,./logs/atom_memcheck.log) $(call escape,$(CURDIR)/atom_test$(EXEC) > $(CURDIR)/logs/atom_test.log)
	$(call TESTCOMM,./logs/matrix_memcheck.log) $(call escape,$(CURDIR)/matrix_test$(EXEC) > $(CURDIR)/logs/matrix_test.log)
	$(call TESTCOMM,./logs/matrix_arithmetic_memcheck.log) $(call escape,$(CURDIR)/matrix_arithmetic_test$(EXEC) > $(CURDIR)/logs/matrix_arithmetic_test.log)
	$(call TESTCOMM,./logs/eigenvalue_memcheck.log) $(call escape,$(CURDIR)/eigenvalue_test$(EXEC) > $(CURDIR)/logs/eigenvalue_test.log)
	$(call TESTCOMM,./logs/geometry_memcheck.log) $(call escape,$(CURDIR)/geometry_test$(EXEC) > $(CURDIR)/logs/geometry_test.log)
	$(call TESTCOMM,./logs/hessian_memcheck.log) $(call escape,$(CURDIR)/hessian_test$(EXEC) > $(CURDIR)/logs/hessian_test.log)
	
.PHONY: clean
clean:
	$(RM) $(EXECS)
ifneq ($(wildcard $(CURDIR)/logs/*.log),)
	$(RM) $(call escape,$(wildcard $(CURDIR)/logs/*.log))
endif


test_test: test_test.cpp test.hpp

atom_test: atom_test.cpp ../Molecule/molecule.hpp test.hpp

matrix_test: matrix_test.cpp ../Base/matrix.hpp ../Base/matrix.cpp

matrix_arithmetic_test: matrix_arithmetic_test.cpp ../Base/matrix_arithmetic_default.hpp ../Base/matrix_arithmetic_default.cpp ../Base/math.o

eigenvalue_test: eigenvalue_test.cpp ../Base/matrix.hpp ../Base/math.hpp ../Base/eigenvalues.hpp \
	../Base/eigenvalues_default.hpp ../Base/math.cpp ../Base/eigenvalues_default_double.cpp ../Base/matrix.cpp

geometry_test: geometry_test.cpp ../Project-1/geometry_default.hpp \
	test.hpp ../Molecule/molecule.hpp ../Molecule/molecule.cpp ../Project-1/geometry_default.cpp \
	../Molecule/molecule_default.hpp ../Base/matrix.hpp ../Base/eigenvalues.hpp ../Base/math.cpp \
	../Base/eigenvalues_default.hpp ../Base/eigenvalues_default_double.cpp ../Base/matrix.cpp
	$(CXX) $(CXXFLAGS) -o geometry_test $^ $(LDFLAGS) $(LDLIBS)
	
hessian_test: hessian_test.cpp ../Project-2/hessian_default.hpp test.hpp ../Molecule/molecule_default.hpp \
	../Base/matrix_default.hpp ../Base/eigenvalues_default.hpp ../Base/matrix.cpp ../Base/eigenvalues_default_double.cpp \
	../Base/math.cpp ../Molecule/molecule.cpp

scf_test: scf_test.cpp ../Project-3/scf_default.hpp ../Project-3/scf_default.cpp test.hpp ../Molecule/molecule_default.hpp \
	../Base/matrix_default.hpp ../Base/eigenvalues_default.hpp ../Base/matrix.cpp ../Base/eigenvalues_default_double.cpp \
	../Base/math.cpp ../Molecule/molecule.cpp ../Base/matrix_arithmetic_default.hpp ../Base/matrix_arithmetic_default.cpp