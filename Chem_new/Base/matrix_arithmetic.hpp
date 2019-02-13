/*
 * matrix_arithmetic.hpp
 *
 *  Created on: Feb 4, 2019
 *      Author: Connor
 */

#ifndef BASE_MATRIX_ARITHMETIC_HPP_
#define BASE_MATRIX_ARITHMETIC_HPP_

#include "matrix.hpp"

namespace compchem {

template<typename _T>
class MatrixArithmeticStrategy {
public:
	MatrixArithmeticStrategy() {
		;
	}
	virtual ~MatrixArithmeticStrategy() = default;

	virtual compchem::AbstractMatrix<_T> &add(
	        const compchem::AbstractMatrix<_T> &a,
	        const compchem::AbstractMatrix<_T> &b) = 0;
	virtual compchem::AbstractMatrix<_T> &subtract(
	        const compchem::AbstractMatrix<_T> &a,
	        const compchem::AbstractMatrix<_T> &b) = 0;
	virtual compchem::AbstractMatrix<_T> &mult(
	        const compchem::AbstractMatrix<_T> &a,
	        const compchem::AbstractMatrix<_T> &b) = 0;
	virtual compchem::AbstractMatrix<_T> &mult(
			const compchem::AbstractMatrix<_T> &a,
			const _T &b) = 0;
	virtual compchem::AbstractMatrix<_T> &inverse(
	        const compchem::AbstractMatrix<_T> &a) = 0;
	virtual compchem::AbstractMatrix<_T> &transpose(
	        const compchem::AbstractMatrix<_T> &a) = 0;

	//Returns the rank, and puts the resulting constants in the consts array,
	//and the resulting coefficients in the coeffs array.
	virtual int solve(compchem::AbstractMatrix<_T> &coeffs,
	        compchem::AbstractMatrix<_T> &consts) = 0;

	virtual _T det(const compchem::AbstractMatrix<_T> &mat) = 0;
	virtual _T trace(const compchem::AbstractMatrix<_T> &mat) = 0;
};

}

#endif /* BASE_MATRIX_ARITHMETIC_HPP_ */
