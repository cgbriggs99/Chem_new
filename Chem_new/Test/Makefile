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
	eigenvalue_test$(EXEC) geometry_test$(EXEC) hessian_test$(EXEC) scf_test$(EXEC) mp2_test$(EXEC) \
	ccsd_test$(EXEC) ccsdt_test$(EXEC) diis_scf_test$(EXEC)

.PHONY: all
all: $(EXECS)
	$(info $^)

%$(EXEC): %.cpp
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -o $@ $^ $(LDLIBS)

#Be aware, some of these tests (scf, mp2) can take a long time. Run at your own risk.
.PHONY: run_tests
run_tests: $(EXECS)
	$(call TESTCOMM,./logs/test_memcheck.log) $(call escape,$(CURDIR)/test_test$(EXEC) > $(CURDIR)/logs/test_test.log)
	$(call TESTCOMM,./logs/atom_memcheck.log) $(call escape,$(CURDIR)/atom_test$(EXEC) > $(CURDIR)/logs/atom_test.log)
	$(call TESTCOMM,./logs/matrix_memcheck.log) $(call escape,$(CURDIR)/matrix_test$(EXEC) > $(CURDIR)/logs/matrix_test.log)
	$(call TESTCOMM,./logs/matrix_arithmetic_memcheck.log) $(call escape,$(CURDIR)/matrix_arithmetic_test$(EXEC) > $(CURDIR)/logs/matrix_arithmetic_test.log)
	$(call TESTCOMM,./logs/eigenvalue_memcheck.log) $(call escape,$(CURDIR)/eigenvalue_test$(EXEC) > $(CURDIR)/logs/eigenvalue_test.log)
	$(call TESTCOMM,./logs/geometry_memcheck.log) $(call escape,$(CURDIR)/geometry_test$(EXEC) > $(CURDIR)/logs/geometry_test.log)
	$(call TESTCOMM,./logs/hessian_memcheck.log) $(call escape,$(CURDIR)/hessian_test$(EXEC) > $(CURDIR)/logs/hessian_test.log)
#	$(call TESTCOMM,./logs/scf_memcheck.log) $(call escape,$(CURDIR)/scf_test$(EXEC) > $(CURDIR)/logs/scf_test.log)
#	$(call TESTCOMM,./logs/mp2_memcheck.log) $(call escape,$(CURDIR)/mp2_test$(EXEC) > $(CURDIR)/logs/mp2_test.log)
#	$(call TESTCOMM,./logs/ccsd_memcheck.log) $(call escape,$(CURDIR)/ccsd_test$(EXEC) > $(CURDIR)/logs/ccsd_test.log)
#	$(call TESTCOMM,./logs/ccsdt_memcheck.log) $(call escape,$(CURDIR)/ccsdt_test$(EXEC) > $(CURDIR)/logs/ccsdt_test.log)
	$(call TESTCOMM,./logs/diis_scf_memcheck.log) $(call escape,$(CURDIR)/diis_scf_test$(EXEC) > $(CURDIR)/logs/diis_scf_test.log)
	
.PHONY: clean
clean:
	$(RM) $(EXECS)
ifneq ($(wildcard $(CURDIR)/logs/*.log),)
	$(RM) $(call escape,$(wildcard $(CURDIR)/logs/*.log))
endif


test_test$(EXEC): test_test.cpp test.hpp

atom_test$(EXEC): atom_test.cpp ../Molecule/molecule.hpp test.hpp

matrix_test$(EXEC): matrix_test.cpp ../Base/matrix.hpp ../Base/matrix.cpp

matrix_arithmetic_test$(EXEC): matrix_arithmetic_test.cpp ../Base/matrix_arithmetic_default.hpp ../Base/matrix_arithmetic_default.cpp ../Base/math.o

eigenvalue_test$(EXEC): eigenvalue_test.cpp ../Base/matrix.hpp ../Base/math.hpp ../Base/eigenvalues.hpp \
	../Base/eigenvalues_default.hpp ../Base/math.cpp ../Base/eigenvalues_default_double.cpp ../Base/matrix.cpp

geometry_test$(EXEC): geometry_test.cpp ../Project-1/geometry_default.hpp \
	test.hpp ../Molecule/molecule.hpp ../Molecule/molecule.cpp ../Project-1/geometry_default.cpp \
	../Molecule/molecule_default.hpp ../Base/matrix.hpp ../Base/eigenvalues.hpp ../Base/math.cpp \
	../Base/eigenvalues_default.hpp ../Base/eigenvalues_default_double.cpp ../Base/matrix.cpp

hessian_test$(EXEC): hessian_test.cpp ../Project-2/hessian_default.hpp test.hpp ../Molecule/molecule_default.hpp \
	../Base/matrix_default.hpp ../Base/eigenvalues_default.hpp ../Base/matrix.cpp ../Base/eigenvalues_default_double.cpp \
	../Base/math.cpp ../Molecule/molecule.cpp

scf_test$(EXEC): scf_test.cpp ../Project-3/scf_default.hpp ../Project-3/scf_default.cpp test.hpp ../Molecule/molecule_default.hpp \
	../Base/matrix_default.hpp ../Base/eigenvalues_default.hpp ../Base/matrix.cpp ../Base/eigenvalues_default_double.cpp \
	../Base/math.cpp ../Molecule/molecule.cpp ../Base/matrix_arithmetic_default.hpp ../Base/matrix_arithmetic_default.cpp

mp2_test$(EXEC): mp2_test.cpp ../Project-4/mp2_default.cpp ../Project-4/mp2_default.hpp ../Project-3/scf_default.hpp \
	../Project-3/scf_default.cpp ../Molecule/molecule_default.hpp ../Base/matrix_default.hpp ../Base/eigenvalues_default.hpp\
	../Base/matrix.cpp ../Base/eigenvalues_default_double.cpp ../Base/math.cpp ../Molecule/molecule.cpp \
	../Base/matrix_arithmetic_default.hpp ../Base/matrix_arithmetic_default.cpp test.hpp

ccsd_test$(EXEC): ccsd_test.cpp ../Project-5/ccsd_default.cpp ../Project-5/ccsd_default.hpp ../Project-3/scf_default.hpp \
	../Project-3/scf_default.cpp ../Molecule/molecule_default.hpp ../Base/matrix_default.hpp ../Base/eigenvalues_default.hpp\
	../Base/matrix.cpp ../Base/eigenvalues_default_double.cpp ../Base/math.cpp ../Molecule/molecule.cpp \
	../Base/matrix_arithmetic_default.hpp ../Base/matrix_arithmetic_default.cpp test.hpp

ccsdt_test$(EXEC): ccsdt_test.cpp ../Project-6/ccsdt.cpp ../Project-6/ccsdt.hpp ../Project-3/scf_default.hpp \
	../Project-3/scf_default.cpp ../Molecule/molecule_default.hpp ../Base/matrix_default.hpp ../Base/eigenvalues_default.hpp\
	../Base/matrix.cpp ../Base/eigenvalues_default_double.cpp ../Base/math.cpp ../Molecule/molecule.cpp \
	../Base/matrix_arithmetic_default.hpp ../Base/matrix_arithmetic_default.cpp test.hpp

diis_scf_test$(EXEC): diis_scf_test.cpp ../Project-8/diis_scf.hpp ../Project-8/diis_scf_default.cpp \
	../Project-3/scf_default.hpp ../Project-3/scf_default.cpp test.hpp ../Molecule/molecule_default.hpp \
	../Base/matrix_default.hpp ../Base/eigenvalues_default.hpp ../Base/matrix.cpp ../Base/eigenvalues_default_double.cpp \
	../Base/math.cpp ../Molecule/molecule.cpp ../Base/matrix_arithmetic_default.hpp ../Base/matrix_arithmetic_default.cpp \
	../Base/matrix_tei.hpp ../Base/matrix_tei.cpp