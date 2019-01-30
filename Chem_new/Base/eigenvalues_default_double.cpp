/*
 * eigenvalues_default_double.cpp
 *
 *  Created on: Jan 28, 2019
 *      Author: connor
 */

#include "eigenvalues_default.hpp"
#include "../Base/matrix.hpp"

#include <stdlib.h>
#include <exception>
#include <stdexcept>

template<>
compchem::Matrix<double> &compchem::strategies::LapackEigenvalues<double>::eigenvals(
		const compchem::Matrix<double> &mat) {
	double *work = (double *) malloc((size_t) (mat.getSize() * sizeof(double)));
	double *out = (double *) malloc(mat.getShape(0) * sizeof(double));
	double *temp = (double *) malloc(mat.getShape(0) * sizeof(double));
	for (int i = 0; i < mat.getShape(0); i++) {
		for (int j = 0; j < mat.getShape(1); j++) {
			work[i * mat.getShape(1) + j] = mat.getEntry( { i, j });
		}
	}
	int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'N', mat.getShape(0), work,
			mat.getShape(1), out, temp, nullptr, mat.getShape(0), nullptr,
			mat.getShape(1));
	free(work);
	Matrix<double> &retval = *new Matrix<double>(out, { mat.getShape(0) });
	free(out);
	free(temp);

	if (info != 0) {
		printf("Error %d\n", info);
		throw(new std::runtime_error("Lapacke error!"));
	}
	return (retval);
}

template<>
compchem::Matrix<double> &compchem::strategies::LapackEigenvalues<double>::eigenvecs_left(
		const compchem::Matrix<double> &mat) {
	double *work = (double *) malloc((size_t) (mat.getSize() * sizeof(double)));
	double *temp1 = new double[mat.getShape(0)];
	double *temp2 = new double[mat.getShape(0)];
	double *out = new double[mat.getSize()];
	for (int i = 0; i < mat.getShape(0); i++) {
		for (int j = 0; j < mat.getShape(1); j++) {
			work[i * mat.getShape(1) + j] = mat.getEntry( { i, j });
		}
	}
	LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'V', 'N', mat.getShape(0), work,
			mat.getShape(1), temp1, temp2, out, mat.getShape(1), nullptr,
			mat.getShape(1));
	free(work);
	Matrix<double> &retval = *new Matrix<double>(out,
			{ mat.getShape(0), mat.getShape(1) });
	delete out;
	delete temp1;
	delete temp2;
	return (retval);
}

template<>
compchem::Matrix<double> &compchem::strategies::LapackEigenvalues<double>::eigenvecs_right(
		const compchem::Matrix<double> &mat) {
	double *work = (double *) malloc((size_t) (mat.getSize() * sizeof(double)));
	double *temp1 = new double[mat.getShape(0)];
	double *temp2 = new double[mat.getShape(0)];
	double *out = new double[mat.getSize()];
	for (int i = 0; i < mat.getShape(0); i++) {
		for (int j = 0; j < mat.getShape(1); j++) {
			work[i * mat.getShape(1) + j] = mat.getEntry( { i, j });
		}
	}
	LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'V', mat.getShape(0), work,
			mat.getShape(1), temp1, temp2, nullptr, mat.getShape(0), out,
			mat.getShape(1));
	free(work);
	Matrix<double> &retval = *new Matrix<double>(out,
			{ mat.getShape(0), mat.getShape(1) });
	delete out;
	delete temp1;
	delete temp2;
	return (retval);
}

template<>
void compchem::strategies::LapackEigenvalues<double>::eigen_all(
		const compchem::Matrix<double> &mat, compchem::Matrix<double> *&evals,
		compchem::Matrix<double> *&rvecs, compchem::Matrix<double> *&lvecs) {
	double *work = (double *) malloc((size_t) (mat.getSize() * sizeof(double)));

	double *outvals = new double[mat.getShape(0)];
	double *temp = new double[mat.getShape(0)];
	double *outvl = new double[mat.getSize()];
	double *outvr = new double[mat.getSize()];
	for (int i = 0; i < mat.getShape(0); i++) {
		for (int j = 0; j < mat.getShape(1); j++) {
			work[i * mat.getShape(1) + j] = mat.getEntry( { i, j });
		}
	}
	LAPACKE_dgeev(LAPACK_ROW_MAJOR, 'N', 'N', mat.getShape(0), work,
			mat.getShape(1), outvals, temp, outvl, mat.getShape(1), outvr,
			mat.getShape(1));
	free(work);

	*evals = *new Matrix<double>(outvals, { mat.getShape(0) });
	*rvecs = *new Matrix<double>(outvr, { mat.getShape(0), mat.getShape(1) });
	*lvecs = *new Matrix<double>(outvl, { mat.getShape(0), mat.getShape(1) });
	delete outvals;
	delete temp;
	delete rvecs;
	delete lvecs;
}

