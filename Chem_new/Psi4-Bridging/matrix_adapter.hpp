/*
 * matrix_bridge.hpp
 *
 *  Created on: Feb 20, 2019
 *      Author: connor
 */

#ifndef MATRIX_ADAPTER_HPP_
#define MATRIX_ADAPTER_HPP_

#include "../Base/base.hpp"

namespace psi {
class Matrix;
typedef std::shared_ptr<Matrix> SharedMatrix;
class Vector3;
}

namespace compchem {

class Psi4MatrixAdapter {
protected:
	Psi4MatrixAdapter() {
		;
	}
	static Psi4MatrixAdapter *singleton;

public:
	virtual ~Psi4MatrixAdapter() {
		;
	}

	static const Psi4MatrixAdapter *getSingleton() {
		if(singleton == nullptr) {
			singleton = new Psi4MatrixAdapter();
		}
		return (singleton);
	}

	static void init() {
		singleton = new Psi4MatrixAdapter();
	}

	virtual compchem::Matrix<double> &convert(psi::Matrix &mat) const;
	virtual compchem::Matrix<double> &convert(psi::SharedMatrix mat) const;
	virtual compchem::Matrix<double> &convert(psi::Vector3 &vec) const;
};

}

#include "matrix_adapter.cpp"

#endif /* MATRIX_ADAPTER_HPP_ */
