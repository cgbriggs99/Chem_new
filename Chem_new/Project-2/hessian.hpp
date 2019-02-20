/*
 * hessian.hpp
 *
 *  Created on: Feb 4, 2019
 *      Author: connor
 */

#ifndef PROJECT_2_HESSIAN_HPP_
#define PROJECT_2_HESSIAN_HPP_

#include <vector>
#include "../Base/matrix.hpp"
#include "../Molecule/molecule.hpp"

namespace compchem {

class HessianStrategy {
public:
	HessianStrategy() {
		;
	}
	virtual ~HessianStrategy() = default;

	virtual compchem::AbstractMatrix<double> &massWeightMatrix(const compchem::AbstractMolecule &mol, const compchem::AbstractMatrix<double> &hess) = 0;
	virtual compchem::AbstractMatrix<double> &computeHessianFreqs(const compchem::AbstractMatrix<double> &hessian_eigs) = 0;
};

}

#endif /* PROJECT_2_HESSIAN_HPP_ */
