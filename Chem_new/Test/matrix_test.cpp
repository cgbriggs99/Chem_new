/*
 * matrix_test.cpp
 *
 *  Created on: Dec 22, 2018
 *      Author: connor
 */

#include "../Base/matrix_default.hpp"
#include "test.hpp"
#include <time.h>
#include <stdlib.h>
#include "../Project-11/disk_allocator.hpp"

class EchoAllocation {
public:
	EchoAllocation() {
		puts("Allocated!");
	}

	~EchoAllocation() {
		puts("Deallocated!");
	}
};

class MatrixTest : public test::Test {
private:
	compchem::Matrix<int> *mat;
	compchem::Matrix<int, compchem::strategies::DiskAllocator<int> > *disk;
public:
	MatrixTest() {
		mat = nullptr;
		disk = nullptr;
		clock_t seed = clock();
		srand(seed);
	}
	~MatrixTest() {
		;
	}
	void testAllocation() {
		mat = new compchem::Matrix<int>( {100, 100, 100});
		try {
			assert(mat->getSize() == 100 * 100 * 100);
		} catch(test::AssertionFailedException *e) {
			printf("Failed! %d successes, %d failures, matrix size is %d.\n",
			        e->getSuccesses(), e->getFails(), mat->getSize());
		}
	}

	void testDiskAllocation() {
		disk = new compchem::Matrix<int,
		        compchem::strategies::DiskAllocator<int> >( {100, 100, 100});
		assert(disk->getSize() == 100 * 100 * 100);
	}

	void testDimensions() {
		assert(mat->getDimension() == 3);
		assert(disk->getDimension() == 3);
	}

	void testStorageRecall() {
		for(int i = 0; i < 100 * 100 * 100; i++) {
			try {
				mat->setEntry(i, i / 10000, (i % 10000) / 100, i % 100);
			} catch(std::out_of_range *e) {
				puts(e->what());
				exit(-1);
			}
		}

		for(int i = 0; i < 1000; i++) {
			int index = (int) rand();
			if(index >= mat->getSize()) {
				index %= mat->getSize();
			}
			try {
				assert(
				        mat->getEntry(
				                index / 10000, (index % 10000) / 100, index
				                        % 100) == index);
			} catch(std::out_of_range *e) {
				puts(e->what());
				exit(-1);
			} catch(test::AssertionFailedException *e) {
				printf("Assertion failed! Got %d at entry %d.\n",
				        mat->getEntry(
				                index / 10000, (index % 10000) / 100, index
				                        % 100), index);
				exit(-1);
			}
		}
	}

	void testDiskStorageRecall() {
		for(int i = 0; i < 100 * 100 * 100; i++) {
			try {
				disk->setEntry(i, i / 10000, (i % 10000) / 100, i % 100);
			} catch(std::out_of_range *e) {
				puts(e->what());
				exit(-1);
			}
		}

		for(int i = 0; i < 1000; i++) {
			int index = (int) rand();
			if(index >= disk->getSize()) {
				index %= disk->getSize();
			}
			try {
				assert(
				        disk->getEntry(
				                index / 10000, (index % 10000) / 100, index
				                        % 100) == index);
			} catch(std::out_of_range *e) {
				puts(e->what());
				exit(-1);
			} catch(test::AssertionFailedException *e) {
				printf("Assertion failed! Got %d at entry %d.\n",
				        disk->getEntry(
				                index / 10000, (index % 10000) / 100, index
				                        % 100), index);
				exit(-1);
			}
		}
	}

	void testDestruction() {
		delete mat;
		delete disk;
	}

	void runTest() override {
		testAllocation();
		testDiskAllocation();
		testDimensions();
		testStorageRecall();
		testDiskStorageRecall();
		testDestruction();
	}
};

int main(void) {
	MatrixTest test;
	test.runTest();
	return (0);
}
