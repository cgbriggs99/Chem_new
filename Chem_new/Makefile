include delete.mk

.PHONY: all
all: link_dlls molecule test base project1

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
	
.PHONY: link_dlls
link_dlls:
	$(MAKE) -C ./Windows all
	
.PHONY: clean
clean:
	$(MAKE) -C ./Test clean
	$(MAKE) -C ./Base clean
	$(MAKE) -C ./Molecule clean
	$(MAKE) -C ./Project-1 clean
	$(MAKE) -C ./Windows clean