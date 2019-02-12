/*
 * matrix_arithmetic_default_double.cpp
 *
 *  Created on: Feb 4, 2019
 *      Author: Connor
 */

#include "matrix_arithmetic_default.hpp"
#include "../Base/matrix_default.hpp"

#ifndef __MATRIX_ARITHMETIC_DEFAULT_CPP__
#define __MATRIX_ARITHMETIC_DEFAULT_CPP__

template<typename _T>
compchem::AbstractMatrix<_T> &compchem::strategies::DefaultMatrixArithmeticStrategy<_T>::add(
		const compchem::AbstractMatrix<_T> &a, const compchem::AbstractMatrix<_T> &b) {
	compchem::Matrix<_T> *out = new compchem::Matrix<_T>(a);
	for(int i = 0; i < out->getShape(0); i++) {
		for(int j = 0; j < out->getShape(1); j++) {
			out->setEntry(a.getEntry(i, j) + b.getEntry(i, j), i, j);
		}
	}
	return (*out);
}

template<typename _T>
compchem::AbstractMatrix<_T> &compchem::strategies::DefaultMatrixArithmeticStrategy<_T>::subtract(const compchem::AbstractMatrix<_T> &a, const compchem::AbstractMatrix<_T> &b) {
	compchem::Matrix<_T> *out = new compchem::Matrix<_T>(a);
	for(int i = 0; i < out->getShape(0); i++) {
		for(int j = 0; j < out->getShape(1); j++) {
			out->setEntry(a.getEntry(i, j) - b.getEntry(i, j), i, j);
		}
	}
	return (*out);
}

template<typename _T>
compchem::AbstractMatrix<_T> &compchem::strategies::DefaultMatrixArithmeticStrategy<_T>::mult(
		const compchem::AbstractMatrix<_T> &a, const compchem::AbstractMatrix<_T> &b) {
	compchem::Matrix<_T> *out = new compchem::Matrix<_T>({a.getShape(0), b.getShape(1)});

	for(int i = 0; i < a.getShape(0); i++) {
		for(int j = 0; j < b.getShape(1); j++) {
			double sum = 0;
			for(int k = 0; k < b.getShape(0); k++) {
				sum += a.getEntry(i, k) * b.getEntry(k, j);
			}
			out->setEntry(sum, i, j);
		}
	}
	return (*out);
}

#define abs(x) ((x < 0)? -x: x)

template<typename _T>
compchem::AbstractMatrix<_T> &compchem::strategies::DefaultMatrixArithmeticStrategy<_T>::inverse(const compchem::AbstractMatrix<_T> &a) {
	compchem::Matrix<_T> *out = new compchem::Matrix<_T>({a.getShape(0), a.getShape(0)}), *work = new compchem::Matrix<_T>(a);

	for(int i = 0; i < out->getShape(0); i++) {
		for(int j = 0; j < out->getShape(1); j++) {
			out->setEntry(0, i, j);
		}
		out->setEntry(1, i, i);
	}

	//Gaussian elimination.
	for(int i = 0; i < work->getShape(0); i++) {
		//Search for absolute largest element.
		_T max = 0;
		int max_ind = 0;
		for(int j = i; j < work->getShape(1); j++) {
			if(abs(work->getEntry(j, i)) > abs(max)) {
				max_ind = j;
				max = work->getEntry(j, i);
			}
		}

		//Swap to put max into the diagonal.
		for(int j = 0; j < work->getShape(1); j++) {
			_T swap1 = work->getEntry(i, j);
			_T swap2 = out->getEntry(i, j);
			work->setEntry(work->getEntry(max_ind, j), i, j);
			out->setEntry(out->getEntry(max_ind, j), i, j);
			work->setEntry(swap1, max_ind, j);
			out->setEntry(swap2, max_ind, j);
		}

		//Divide the row with the max.
		for(int j = 0; j < work->getShape(1); j++) {
			work->setEntry(work->getEntry(i, j) / max, i, j);
			out->setEntry(out->getEntry(i, j) / max, i, j);
		}

		//Cancel all the others in this row.
		for(int j = 0; j < work->getShape(0); j++) {
			if(j == i) {
				continue;
			}
			_T coef = work->getEntry(j, i);
			for(int k = 0; k < work->getShape(1); k++) {
				work->setEntry(work->getEntry(j, k) - coef * work->getEntry(i, k), j, k);
				out->setEntry(out->getEntry(j, k) - coef * out->getEntry(i, k), j, k);
			}
		}
	}
	delete work;
	return (*out);
}

template<typename _T>
compchem::AbstractMatrix<_T> &compchem::strategies::DefaultMatrixArithmeticStrategy<_T>::transpose(const compchem::AbstractMatrix<_T> &a) {
	compchem::Matrix<_T> *out = new compchem::Matrix<_T>({a.getShape(1), a.getShape(0)});
	for(int i = 0; i < a.getShape(0); i++) {
		for(int j = 0; j < a.getShape(1); j++) {
			out->setEntry(a.getEntry(i, j), j, i);
		}
	}
	return (*out);
}


#endif
