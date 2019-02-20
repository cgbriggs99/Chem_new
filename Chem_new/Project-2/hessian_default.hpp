/*
 * hessian_default.hpp
 *
 *  Created on: Feb 4, 2019
 *      Author: connor
 */

#ifndef PROJECT_2_HESSIAN_DEFAULT_HPP_
#define PROJECT_2_HESSIAN_DEFAULT_HPP_

#include "hessian.hpp"
#include "../Base/matrix_default.hpp"
#include <vector>
#include <cmath>

namespace compchem {
namespace strategies {
class DefaultHessianStrategy : public HessianStrategy {
public:
	DefaultHessianStrategy() {
		;
	}
	~DefaultHessianStrategy() {
		;
	}

	compchem::AbstractMatrix<double> &massWeightMatrix(const compchem::AbstractMolecule &mol,
			const compchem::AbstractMatrix<double> &hess) {
		compchem::Matrix<double> *out = new compchem::Matrix<double>({mol.natom() * 3, mol.natom() * 3});
		for(int i = 0; i < mol.natom() * 3; i++) {
			for(int j = 0; j < mol.natom() * 3; j++) {
				out->setEntry(hess.getEntry(i, j) / sqrt(mol.mass(i / 3) * mol.mass(j / 3)), i, j);
			}
		}
		return (*out);
	}

	compchem::AbstractMatrix<double> &computeHessianFreqs(const compchem::AbstractMatrix<double> &hessian_eigs) override {
		compchem::Matrix<double> *out = new compchem::Matrix<double>({hessian_eigs.getSize()});
		for(int i = 0; i < hessian_eigs.getSize(); i++) {
			double hold = hessian_eigs.getEntry(i);
			if(hold <= 0 || hold == -0.0) {
				out->setEntry(0, i);
			} else {
				out->setEntry(sqrt(hessian_eigs.getEntry(i)) * 5140.697669352, i);
			}
		}
		return (*out);
	}
};

}
}

#endif /* PROJECT_2_HESSIAN_DEFAULT_HPP_ */
