/*
 * matrix_arithmetic_default.hpp
 *
 *  Created on: Feb 4, 2019
 *      Author: Connor
 */

#ifndef BASE_MATRIX_ARITHMETIC_DEFAULT_HPP_
#define BASE_MATRIX_ARITHMETIC_DEFAULT_HPP_

#include "matrix_arithmetic.hpp"


namespace compchem {
namespace strategies {

template<typename _T>
class DefaultMatrixArithmeticStrategy : public compchem::MatrixArithmeticStrategy<_T> {
public:
	DefaultMatrixArithmeticStrategy() {
		;
	}

	virtual ~DefaultMatrixArithmeticStrategy() {
		;
	}

	compchem::AbstractMatrix<_T> &add(const compchem::AbstractMatrix<_T> &a, const compchem::AbstractMatrix<_T> &b) override;
	compchem::AbstractMatrix<_T> &subtract(const compchem::AbstractMatrix<_T> &a, const compchem::AbstractMatrix<_T> &b) override;
	compchem::AbstractMatrix<_T> &mult(const compchem::AbstractMatrix<_T> &a, const compchem::AbstractMatrix<_T> &b) override;
	compchem::AbstractMatrix<_T> &inverse(const compchem::AbstractMatrix<_T> &a) override;
	compchem::AbstractMatrix<_T> &transpose(const compchem::AbstractMatrix<_T> &a) override;
	int solve(compchem::AbstractMatrix<_T> &coeffs,
				compchem::AbstractMatrix<_T> &consts) override;

	_T det(const compchem::AbstractMatrix<_T> &mat) override;
	_T trace(const compchem::AbstractMatrix<_T> &mat) override;
};

}
}

#include "matrix_arithmetic_default.cpp"

#endif /* BASE_MATRIX_ARITHMETIC_DEFAULT_HPP_ */
