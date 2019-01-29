

.PHONY: all
all: molecule test base project1

.PHONY: test
test: 
	$(MAKE) -C ./Test all 
	
.PHONY: molecule
molecule:
	$(MAKE) -C ./Molecule all
	
.PHONY: base
base:
	$(MAKE) -C ./Base all
	
.PHONY: project1
project1:
	$(MAKE) -C ./Project-1 all

.PHONY: run_tests
run_tests:
	$(MAKE) -C ./Test run_tests
	
.PHONY: clean
clean:
	rm -rf *.o
	$(MAKE) -C ./Test clean
	$(MAKE) -C ./Base clean
	$(MAKE) -C ./Molecule clean
	$(MAKE) -C ./Project-1 clean