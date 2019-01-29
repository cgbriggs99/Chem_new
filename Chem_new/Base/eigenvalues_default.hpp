/*
 * eigenvalues_default.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: connor
 */

#ifndef BASE_EIGENVALUES_DEFAULT_HPP_
#define BASE_EIGENVALUES_DEFAULT_HPP_

#include <complex>
#include "eigenvalues.hpp"
#include <lapacke.h>

namespace compchem {
namespace strategies {


template<typename T>
class LapackEigenvalues : public compchem::EigenvalueStrategy<T> {
public:
	LapackEigenvalues() {
		;
	}

	~LapackEigenvalues() {
		;
	}

	Matrix<T> &eigenvals(const Matrix<T> &mat) override;
	Matrix<T> &eigenvecs_left(const Matrix<T> &mat) override;
	Matrix<T> &eigenvecs_right(const Matrix<T> &mat) override;
	void eigen_all(const Matrix<T> &mat, Matrix<T> *&evals, Matrix<T> *&rvecs, Matrix<T> *&lvecs) override;
};




//template<typename T>
//Matrix<T> &eigenvals(const Matrix<T> &mat) {
//	return (NULL);
//}
//
//template<typename T>
//Matrix<T> &eigenvecs_left(const Matrix<T> &mat) {
//	return (NULL);
//}
//
//template<typename T>
//Matrix<T> &eigenvecs_right(const Matrix<T> &mat) {
//	return (NULL);
//}
//
//template<typename T>
//void eigen_all(const Matrix<T> &mat, Matrix<T> *&evals, Matrix<T> *&rvecs, Matrix<T> *&lvecs)  {
//	;
//}

template<>
Matrix<double> &LapackEigenvalues<double>::eigenvals(const Matrix<double> &mat);

template<>
Matrix<double> &LapackEigenvalues<double>::eigenvecs_left(const Matrix<double> &mat);

template<>
Matrix<double> &LapackEigenvalues<double>::eigenvecs_right(const Matrix<double> &mat);

template<>
void LapackEigenvalues<double>::eigen_all(const Matrix<double> &mat,
		Matrix<double> *&evals, Matrix<double> *&rvecs, Matrix<double> *&lvecs);

}
}




#endif /* BASE_EIGENVALUES_DEFAULT_HPP_ */
