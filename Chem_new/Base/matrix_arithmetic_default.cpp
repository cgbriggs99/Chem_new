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
compchem::AbstractMatrix<_T> &compchem::strategies::DefaultMatrixArithmeticStrategy<
        _T>::add(const compchem::AbstractMatrix<_T> &a,
        const compchem::AbstractMatrix<_T> &b) {
	compchem::Matrix<_T> *out = new compchem::Matrix<_T>(a);
	for(int i = 0; i < out->getShape(0); i++) {
		for(int j = 0; j < out->getShape(1); j++) {
			out->setEntry(a.getEntry(i, j) + b.getEntry(i, j), i, j);
		}
	}
	return (*out);
}

template<typename _T>
compchem::AbstractMatrix<_T> &compchem::strategies::DefaultMatrixArithmeticStrategy<
        _T>::subtract(const compchem::AbstractMatrix<_T> &a,
        const compchem::AbstractMatrix<_T> &b) {
	compchem::Matrix<_T> *out = new compchem::Matrix<_T>(a);
	for(int i = 0; i < out->getShape(0); i++) {
		for(int j = 0; j < out->getShape(1); j++) {
			out->setEntry(a.getEntry(i, j) - b.getEntry(i, j), i, j);
		}
	}
	return (*out);
}

template<typename _T>
compchem::AbstractMatrix<_T> &compchem::strategies::DefaultMatrixArithmeticStrategy<
        _T>::mult(const compchem::AbstractMatrix<_T> &a,
        const compchem::AbstractMatrix<_T> &b) {
	if(a.getShape(1) != b.getShape(0)) {
		throw(new std::length_error(
		        "Error! Mismatched dimensions in matrices to be multiplied!"));
	}

	compchem::Matrix<_T> *out = new compchem::Matrix<_T>(
	        {a.getShape(0), b.getShape(1)});

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
compchem::AbstractMatrix<_T> &compchem::strategies::DefaultMatrixArithmeticStrategy<
        _T>::inverse(const compchem::AbstractMatrix<_T> &a) {
	if(a.getShape(0) != a.getShape(1)) {
		throw(std::length_error(
		        "Error! Can only perform inversion on square matrices!"));
	}

	compchem::Matrix<_T> *out = new compchem::Matrix<_T>(
	        {a.getShape(0), a.getShape(0)}), *work = new compchem::Matrix<_T>(
	        a);

	for(int i = 0; i < out->getShape(0); i++) {
		for(int j = 0; j < out->getShape(1); j++) {
			out->setEntry(0, i, j);
		}
		out->setEntry(1, i, i);
	}
	int current_row = 0;

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
			_T swap1 = work->getEntry(current_row, j);
			_T swap2 = out->getEntry(current_row, j);
			work->setEntry(work->getEntry(max_ind, j), current_row, j);
			out->setEntry(out->getEntry(max_ind, j), current_row, j);
			work->setEntry(swap1, max_ind, j);
			out->setEntry(swap2, max_ind, j);
		}

		//Divide the row with the max.
		for(int j = 0; j < work->getShape(1); j++) {
			work->setEntry(work->getEntry(current_row, j) / max, current_row,
			        j);
			out->setEntry(out->getEntry(current_row, j) / max, current_row, j);
		}

		//Cancel all the others in this row.
		for(int j = 0; j < work->getShape(0); j++) {
			if(j == i) {
				continue;
			}
			_T coef = work->getEntry(j, i);
			for(int k = 0; k < work->getShape(1); k++) {
				work->setEntry(
				        work->getEntry(j, k)
				                - coef * work->getEntry(current_row, k), j, k);
				out->setEntry(
				        out->getEntry(j, k)
				                - coef * out->getEntry(current_row, k), j, k);
			}
		}
		current_row++;
	}
	delete work;
	return (*out);
}

template<typename _T>
compchem::AbstractMatrix<_T> &compchem::strategies::DefaultMatrixArithmeticStrategy<
        _T>::transpose(const compchem::AbstractMatrix<_T> &a) {
	compchem::Matrix<_T> *out = new compchem::Matrix<_T>(
	        {a.getShape(1), a.getShape(0)});
	for(int i = 0; i < a.getShape(0); i++) {
		for(int j = 0; j < a.getShape(1); j++) {
			out->setEntry(a.getEntry(i, j), j, i);
		}
	}
	return (*out);
}

template<typename _T>
int compchem::strategies::DefaultMatrixArithmeticStrategy<_T>::solve(
        compchem::AbstractMatrix<_T> &coeffs,
        compchem::AbstractMatrix<_T> &consts) {
	if(coeffs.getShape(0) != consts.getShape(0)) {
		throw(new std::length_error(
		        "Error! Can not perform elimination on matrices with different numbers of rows!"));
	}
	int current_row = 0;

	//Gaussian elimination.
	for(int i = 0; i < std::min(coeffs.getShape(0), coeffs.getShape(1)); i++) {

		//Loop through and find the relative max value.
		_T max = 0;
		int max_ind = current_row;
		for(int j = current_row; j < coeffs.getShape(1); j++) {
			if(abs(coeffs.getEntry(j, i)) > max) {
				max = coeffs.getEntry(j, i);
				max_ind = j;
			}
		}

		//If this row can not be filled,
		if(max == 0) {
			continue;
		}

		if(max_ind != current_row) {
			//Swap the max row to the proper position.
			for(int j = i; j < coeffs.getShape(1); j++) {
				_T swap = coeffs.getEntry(current_row, j);
				coeffs.setEntry(coeffs.getEntry(max_ind, j), current_row, j);
				coeffs.setEntry(swap, max_ind, j);
			}
			for(int j = 0; j < consts.getShape(1); j++) {
				_T swap = consts.getEntry(current_row, j);
				consts.setEntry(consts.getEntry(max_ind, j), current_row, j);
				consts.setEntry(swap, max_ind, j);
			}
		}

		//Scale the max row.
		coeffs.setEntry(1, current_row, i);
		for(int j = i + 1; j < coeffs.getShape(1); j++) {
			coeffs.setEntry(coeffs.getEntry(current_row, j) / max, current_row,
			        j);
		}
		for(int j = 0; j < consts.getShape(1); j++) {
			consts.setEntry(consts.getEntry(current_row, j) / max, current_row,
			        j);
		}

		//Eliminate
		for(int j = 0; j < coeffs.getShape(0); j++) {
			if(j == current_row) {
				continue;
			}
			_T head = coeffs.getEntry(j, i);
			for(int k = i; k < coeffs.getShape(1); k++) {
				coeffs.setEntry(
				        coeffs.getEntry(j, k)
				                - head * coeffs.getEntry(current_row, k), j, k);
			}
			for(int k = 0; k < consts.getShape(1); k++) {
				consts.setEntry(
				        consts.getEntry(j, k)
				                - head * consts.getEntry(current_row, k), j, k);
			}
		}
		current_row++;
	}

	return (current_row);
}

template<typename _T>
_T compchem::strategies::DefaultMatrixArithmeticStrategy<_T>::det(
        const compchem::AbstractMatrix<_T> &mat) {
	if(mat.getShape(0) != mat.getShape(1)) {
		throw(std::length_error(
		        "Error! Can only find the determinant of square matrices!"));
	}

	compchem::Matrix<_T> *work = new compchem::Matrix<_T>(mat);
	_T out = 1;

	int current_row = 0;

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
		if(max_ind != current_row) {
			for(int j = 0; j < work->getShape(1); j++) {
				_T swap1 = work->getEntry(current_row, j);
				work->setEntry(work->getEntry(max_ind, j), current_row, j);
				work->setEntry(swap1, max_ind, j);
			}
			out = -out;
		}

		//Divide the row with the max.
		for(int j = 0; j < work->getShape(1); j++) {
			work->setEntry(work->getEntry(current_row, j) / max, current_row,
			        j);
		}
		out /= max;

		//Cancel all the others in this row.
		for(int j = 0; j < work->getShape(0); j++) {
			if(j == i) {
				continue;
			}
			_T coef = work->getEntry(j, i);
			for(int k = 0; k < work->getShape(1); k++) {
				work->setEntry(
				        work->getEntry(j, k)
				                - coef * work->getEntry(current_row, k), j, k);
			}
		}
		current_row++;
	}
	delete work;
	return (out);
}

template<typename _T>
_T compchem::strategies::DefaultMatrixArithmeticStrategy<_T>::trace(const compchem::AbstractMatrix<_T> &mat) {
	if(mat.getShape(0) != mat.getShape(1)) {
		throw(new std::length_error("Can not find the trace of a non-square matrix."));
	}

	_T sum = 0;
	for(int i = 0; i < mat.getShape(0); i++) {
		sum += mat.getEntry(i, i);
	}
	return (sum);
}

#undef abs

#endif
