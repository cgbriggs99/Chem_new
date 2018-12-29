/*
 * matrix_test.cpp
 *
 *  Created on: Dec 22, 2018
 *      Author: connor
 */

#include "../Base/matrix.hpp"

class TestMatrix {
public:
	static int testInitialization();
	static int testFinalization();
	static int testStorage();
	static int testRetrieval();
	static int testSubclassing();
	static int testAllocators();
};

class SubclassMatrix : public compchem::Matrix<int> {
public:
	SubclassMatrix();
	~SubclassMatrix();

	int &getEntry(std::initializer_list<int> index) override;
	const int &getEntry(std::initializer_list<int> index) const override;
};

class TestAllocator {
	TestAllocator();
	~TestAllocator();

	int *allocate(int n);
	void deallocate(int *ptr, int n);
};
