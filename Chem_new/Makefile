include delete.mk
include psicompile.mk

.PHONY: all
all: link_dlls molecule base project1 test #plugin

.PHONY: test
test: 
	$(info Making the tests)
	$(MAKE) -C ./Test all
	$(info Finished making the tests)
	
.PHONY: rebuild_test
rebuild_test: 
	$(info Making the tests)
	$(MAKE) -B -C ./Test all
	$(info Finished making the tests)
	
.PHONY: molecule
molecule:
	$(info Making molecule things)
	$(MAKE) -C ./Molecule all
	$(info Finished making molecule things)
	
.PHONY: base
base:
	$(info Making utilities.)
	$(MAKE) -C ./Base all
	$(info Finished making utilities)
	
.PHONY: project1
project1:
	$(info Building Project 1.)
	$(MAKE) -C ./Project-1 all
	$(info Finished building Project 1)

.PHONY: run_tests
run_tests:
	$(info Running unit tests.)
	$(MAKE) -C ./Test run_tests
	$(info Finished running unit tests. Output placed in Test/logs directory.)
	
.PHONY: link_dlls
link_dlls:
	$(MAKE) -C ./Windows all
	
.PHONY: plugin
plugin:
	cd ./Plugin
	$(PSI) --plugin-name AOInts
	$(PSI) --plugin-compile
	$(MAKE) install
	cd ../
	
.PHONY: clean
clean:
	$(info Cleaning)
	$(MAKE) -C ./Test clean
	$(MAKE) -C ./Base clean
	$(MAKE) -C ./Molecule clean
	$(MAKE) -C ./Project-1 clean
	$(MAKE) -C ./Windows clean
	$(info Finished cleaning)