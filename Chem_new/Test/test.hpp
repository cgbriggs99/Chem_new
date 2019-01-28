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

class Test {
private:
	int fails;
	int successes;
	int total;
public:
	Test() {
		fails = 0;
		successes = 0;
		total = 0;
	}

	virtual ~Test() {
		if(fails != 0) {
			printf("%d tests failed out of %d with %d successes.\n", fails, total, successes);
		} else {
			printf("%d test executed, all passed.\n", total);
		}
	}

	virtual bool assert(bool cond) {
		total++;
		if(cond) {
			successes++;
		} else {
			fails++;
			throw(new AssertionFailedException(successes, fails));
		}
		return (cond);
	}

	virtual bool assert_print(bool cond) {
		total++;
		try {
			assert(cond);
			printf("Success: s %d f %d\n", successes, fails);
		} catch(AssertionFailedException *e) {
			printf("Failure: s %d f %d\n", successes, fails);
			throw(e);
		}
		return (cond);
	}

	virtual int getSuccesses() {
		return (successes);
	}

	virtual int getFails() {
		return (fails);
	}

	virtual int getTotalTests() {
		return (total);
	}

	virtual void incrementSuccess() {
		successes++;
	}

	virtual void incrementFails() {
		fails++;
	}

	virtual void decrementSuccess() {
		successes--;
	}

	virtual void decrementFails() {
		fails--;
	}

	virtual void runTest() = 0;
};

}

#endif /* TEST_HPP_ */
