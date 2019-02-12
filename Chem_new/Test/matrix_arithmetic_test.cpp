/*
 * matrix_arithmetic_test.cpp
 *
 *  Created on: Feb 4, 2019
 *      Author: Connor
 */

#include "test.hpp"
#include "../Base/matrix_arithmetic_default.hpp"
#include "../Base/math.hpp"

double mat1[] = {
		0.474244, 0.620543, 0.138216, 0.0934553, 0.0265619,
		0.654112, 0.581607, 0.515732, 0.0438419, 0.567722,
		0.793021, 0.890975, 0.107402, 0.522419, 0.218708,
		0.796539, 0.528502, 0.696646, 0.989611, 0.192604,
		0.274031, 0.296648, 0.971876, 0.390136, 0.046910
},
mat1_i[5][5] = {
		6.5641, 1.06937, -6.5381, 4.30869, -3.86695,
		-2.83318, -0.792479, 4.62737, -3.15463, 2.57334,
		0.531452, 0.212138, -0.822414, 0.0101393, 0.924402,
		-3.19146, -0.82274, 2.72224, -0.447714, 0.910555,
		-4.8968, 1.21202, 3.32932, -1.70719, 0.909038
};

class MatrixArithmeticTest : public test::Test {
private:
	compchem::MatrixArithmeticStrategy<double> *strat;

public:
	MatrixArithmeticTest() {
		strat = new compchem::strategies::DefaultMatrixArithmeticStrategy<double>();
	}

	~MatrixArithmeticTest() {
		delete strat;
	}

	void testInverse() {
		compchem::Matrix<double> *mat = new compchem::Matrix<double>(mat1, {5, 5});
		compchem::Matrix<double> *inv = (compchem::Matrix<double> *) &strat->inverse(*mat);
		for(int i = 0; i < 5; i++) {
			for(int j = 0; j < 5; j++) {
				assert(compchem::compareDoubles(inv->getEntry(i, j), mat1_i[i][j], 0.001) == 0);
			}
		}
		delete mat;
		delete inv;
	}

	void runTest() {
		testInverse();
	}
};

int main(void) {
	MatrixArithmeticTest test;

	test.runTest();
	return (0);
}


