/*
 * eigenvalue_test.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: connor
 */

#include "test.hpp"
#include "../Base/eigenvalues.hpp"
#include "../Base/eigenvalues_default.hpp"
#include "../Base/math.hpp"

double __test[] = {
		0.89345, 0.651442, 0.508916, 0.420615, 0.135494,
		0.377386, 0.218746, 0.9852, 0.882663, 0.666169,
		0.13706, 0.651465, 0.951312, 0.547013, 0.150274,
		0.416453, 0.463564, 0.872943, 0.662143, 0.379641,
		0.270043, 0.569471, 0.0576367, 0.983367, 0.239012
};

static int compare(const void *a, const void *b) {
	if((*(double *) a) < (*(double *) b)) {
		return (1);
	} else if((*(double *) a) > (*(double *) b)) {
		return (-1);
	}
	return (0);
}

double values[] = {2.64856, 0.612474, 0.331351, -0.0103405, -0.617385};

class EigenvalueTest : public test::Test {
private:
	compchem::Matrix<double> *mat;
	compchem::strategies::LapackEigenvalues<double> strat;
public:
	EigenvalueTest() {
		mat = new compchem::Matrix<double>(__test, {5, 5});
	}
	~EigenvalueTest() {
		delete mat;
	}

	void testEigenvalues() {
		compchem::Matrix<double> &evecs = strat.eigenvals(*mat);

		double *hold = (double *) malloc(evecs.getSize() * sizeof(double));
		for(int i = 0; i < evecs.getSize(); i++) {
			hold[i] = evecs.getEntry({i});
		}
		qsort(hold, evecs.getSize(), sizeof(double), compare);
		for(int i = 0; i < 5; i++) {
			assert(compchem::compareDoubles(hold[i], values[i], 0.00001) == 0);
		}
		free(hold);
		delete &evecs;
	}

	void runTest() {
		testEigenvalues();
	}
};

int main(void) {
	EigenvalueTest test;
	test.runTest();
	return (0);
}
