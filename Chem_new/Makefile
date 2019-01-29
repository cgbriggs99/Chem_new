

.PHONY: all
all: molecule test

.PHONY: test
test: 
	$(MAKE) -C ./Test all 
	
.PHONY: molecule
molecule:
	$(MAKE) -C ./Molecule all
	
.PHONY: base
base:
	$(MAKE) -C ./Base all

.PHONY: run_tests
run_tests:
	$(MAKE) -C ./Test run_tests
	
.PHONY: clean
clean:
	rm -rf *.o
	$(MAKE) -C ./Test clean