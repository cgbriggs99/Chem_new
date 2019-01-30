/*
 * matrix_test.cpp
 *
 *  Created on: Dec 22, 2018
 *      Author: connor
 */

#include "../Base/matrix.hpp"
#include "test.hpp"
#include <time.h>
#include <stdlib.h>

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
public:
	MatrixTest() {
		mat = nullptr;
		clock_t seed = clock();
		srand(seed);
	}
	~MatrixTest() {
		;
	}
	void testAllocation() {
		mat = new compchem::Matrix<int>({100, 100, 100});
		try {
			assert(mat->getSize() == 100 * 100 * 100);
		} catch(test::AssertionFailedException *e) {
			printf("Failed! %d successes, %d failures, matrix size is %d.\n",
					e->getSuccesses(), e->getFails(), mat->getSize());
		}
	}

	void testDimensions() {
		assert(mat->getDimension() == 3);
	}

	void testStorageRecall() {
		for(int i = 0; i < 100 * 100 * 100; i++) {
			try {
				mat->getEntry({i / 10000, (i % 10000) / 100, i % 100}) = i;
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
				assert(mat->getEntry({index /10000, (index % 10000) / 100, index % 100}) == index);
			} catch(std::out_of_range *e) {
				puts(e->what());
				exit(-1);
			} catch(test::AssertionFailedException *e) {
				printf("Assertion failed! Got %d at entry %d.\n",
						mat->getEntry({index / 10000, (index % 10000) / 100, index % 100}), index);
				exit(-1);
			}
		}
	}

	void testDestruction() {
		delete mat;
	}

	void runTest() override {
		testAllocation();
		testDimensions();
		testStorageRecall();
		testDestruction();
	}
};

int main(void) {
	MatrixTest test;
	test.runTest();
	return (0);
}
