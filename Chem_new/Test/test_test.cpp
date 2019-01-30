/*
 * test_test.cpp
 *
 *  Created on: Jan 28, 2019
 *      Author: connor
 */

#include "test.hpp"

#include <stdio.h>
#include <iostream>

using namespace std;

class TestTest : public test::Test {
public:
	TestTest() :
	Test() {
		;
	}
	~TestTest() {
		;
	}

	void testAssertTrue() {
		puts("Testing assert true...");
		assert(true);
	}

	void testAssertPrintTrue() {
		puts("Testing assert_print true...");
		assert_print(true);
	}

	void testAssertFalse() {
		puts("Testing assert false...");
		try {
			assert(false);
		} catch(test::AssertionFailedException &e) {
			this->incrementSuccess();
			this->decrementFails();
		}
	}

	void testAssertPrintFalse() {
		puts("Testing assert_print false...");
		try {
			assert_print(false);
		} catch(test::AssertionFailedException &e) {
			this->incrementSuccess();
			this->decrementFails();
		}
	}

	void runTest() override {
		testAssertTrue();
		testAssertFalse();
		testAssertPrintTrue();
		testAssertPrintFalse();
	}
};

int main(void) {
	test::Test *test = new TestTest();
	test->runTest();
	delete test;
	return (0);
}

