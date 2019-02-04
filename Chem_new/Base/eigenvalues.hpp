/*
 * eigenvalues.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: connor
 */

#ifndef BASE_EIGENVALUES_HPP_
#define BASE_EIGENVALUES_HPP_

#include "matrix_default.hpp"

namespace compchem {

template<typename T>
class EigenvalueStrategy {
public:
	virtual ~EigenvalueStrategy() = default;

	virtual Matrix<T> &eigenvals(const AbstractMatrix<T> &mat) = 0;
	virtual Matrix<T> &eigenvecs_left(const AbstractMatrix<T> &mat) = 0;
	virtual Matrix<T> &eigenvecs_right(const AbstractMatrix<T> &mat) = 0;
	virtual void eigen_all(const AbstractMatrix<T> &mat, Matrix<T> *&evals, Matrix<T> *&rvecs, Matrix<T> *&lvecs) = 0;
};

}

#endif /* BASE_EIGENVALUES_HPP_ */
