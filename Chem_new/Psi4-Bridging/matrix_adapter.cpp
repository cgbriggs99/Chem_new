/*
 * matrix_adapter.cpp
 *
 *  Created on: Feb 20, 2019
 *      Author: connor
 */

#ifndef __MATRIX_ADAPTER_CPP__
#define __MATRIX_ADAPTER_CPP__

#include "matrix_adapter.hpp"


compchem::Matrix<double> &compchem::Psi4MatrixAdapter::convert(psi::Matrix &mat) const {
	compchem::Matrix<double> *out = new compchem::Matrix<double>({mat.nrow(), mat.ncol()});

	for(int i = 0; i < out->getShape(0); i++) {
		for(int j = 0; j < out->getShape(1); j++) {
			out->setEntry(mat.get(i, j), i, j);
		}
	}
	return (*out);
}

compchem::Matrix<double> &compchem::Psi4MatrixAdapter::convert(psi::SharedMatrix mat) const {
	compchem::Matrix<double> *out = new compchem::Matrix<double>({mat->nrow(), mat->ncol()});

	for(int i = 0; i < out->getShape(0); i++) {
		for(int j = 0; j < out->getShape(1); j++) {
			out->setEntry(mat->get(i, j), i, j);
		}
	}
	return (*out);
}

compchem::Matrix<double> &compchem::Psi4MatrixAdapter::convert(psi::Vector3 &vec) const {
	compchem::Matrix<double> *out = new compchem::Matrix<double>({3});
	out->setEntry(vec.get(0), 0);
	out->setEntry(vec.get(1), 1);
	out->setEntry(vec.get(2), 2);
	return (*out);
}

#endif
