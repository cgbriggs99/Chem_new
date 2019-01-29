/*
 * test.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: connor
 */

#ifndef TEST_HPP_
#define TEST_HPP_

#include <exception>
#include <string>
#include <stdio.h>


namespace test {
//
class AssertionFailedException : public std::exception {
private:
	std::string message;
	int successes, fails;
public:
	AssertionFailedException(int successes, int fails) {
		this->message = "AssertionFailedException!";
		this->successes = successes;
		this->fails = fails;
	}

	const char *what() {
		return (message.c_str());
	}

	int getSuccesses() {
		return (successes);
	}

	int getFails() {
		return (fails);
	}
};

/*
 * Test base class. Unit tests should inherit this class.
 */
class Test {
private:
	int fails;	//number of fails.
	int successes; //number of passes.
	int total;	//number of tests.
public:
	Test() {
		fails = 0;
		successes = 0;
		total = 0;
	}

	/*
	 * Destructor for the test class. It prints the total successes and fails, and if all the tests were successfull.
	 */
	virtual ~Test() {
		if(fails != 0) {
			printf("%d tests failed out of %d with %d successes.\n", fails, total, successes);
		} else {
			printf("%d tests executed, all passed.\n", total);
		}
	}

	/*
	 * This is the basic function for assertions. It takes a condition and, if true, returns and increments
	 * the number of successes and the number of tests, if false, increments the number of fails, throws an
	 * exception and prints.
	 */
	virtual bool assert(bool cond) {
		total++;
		if(cond) {
			successes++;
		} else {
			fails++;
			printf("Failure: s %d f %d\n", successes, fails);
			throw(new AssertionFailedException(successes, fails));
		}
		return (cond);
	}

	/*
	 * This is the basic function for assertions. It takes a condition and, if true, returns and increments
	 * the number of successes and the number of tests, if false, increments the number of fails and throws an
	 * exception.
	 */
	virtual bool assert_noprint(bool cond) {
		total++;
		if(cond) {
			successes++;
		} else {
			fails++;
			throw(new AssertionFailedException(successes, fails));
		}
		return (cond);
	}

	/*
	 * Wrapper for the basic assertion function. This prints the status of the tests, whether it was a success and how
	 * many successes and failures have been met.
	 */
	virtual bool assert_print(bool cond) {
		try {
			assert(cond);
			printf("Success: s %d f %d\n", successes, fails);
		} catch(AssertionFailedException *e) {
			printf("Failure: s %d f %d\n", successes, fails);
			throw(e);
		}
		return (cond);
	}

	/*
	 * Returns the number of successes.
	 */
	virtual int getSuccesses() {
		return (successes);
	}

	/*
	 * Returns the number of failures.
	 */
	virtual int getFails() {
		return (fails);
	}

	/*
	 * Returns the number of tests done so far.
	 */
	virtual int getTotalTests() {
		return (total);
	}

	/*
	 * Increments the number of successes without incrementing the number of tests.
	 */
	virtual void incrementSuccess() {
		successes++;
	}

	/*
	 * Increments the number of fails without incrementing the number of tests.
	 */
	virtual void incrementFails() {
		fails++;
	}

	/*
	 * See incrementSuccess.
	 */
	virtual void decrementSuccess() {
		successes--;
	}

	/*
	 * See incrementFails.
	 */
	virtual void decrementFails() {
		fails--;
	}

	/*
	 * This is the function that should be called to run the tests. It must be overridden in the subclass.
	 */
	virtual void runTest() = 0;
};

}

#endif /* TEST_HPP_ */
