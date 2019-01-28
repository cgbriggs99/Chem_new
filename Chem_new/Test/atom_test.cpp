/*
 * atom_test.cpp
 *
 *  Created on: Jan 28, 2019
 *      Author: connor
 */

#include "../Molecule/molecule.hpp"
#include "test.hpp"

class AtomTest : public test::Test {
public:
	AtomTest() :
		test::Test() {
		;
	}
	~AtomTest() {
		;
	}

	void testStorage() {
		compchem::Atom *hydrogen = new compchem::Atom(1, 1.00789, 0, 0, 0, 0);
		assert(hydrogen->getAtomicNum() == 1);
		assert(hydrogen->getMass() == 1.00789);
		assert(hydrogen->getCharge() == 0);
		assert(hydrogen->getX() == 0);
		assert(hydrogen->getY() == 0);
		assert(hydrogen->getZ() == 0);
	}

	void runTests() override {
		testStorage();
	}
};
