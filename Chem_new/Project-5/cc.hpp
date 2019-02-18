/*
 * ccsd.hpp
 *
 *  Created on: Feb 11, 2019
 *      Author: connor
 */

#ifndef PROJECT_5_CC_HPP_
#define PROJECT_5_CC_HPP_

#include "../Base/base.hpp"

namespace compchem {

class AbstractCCCorrection {
public:
	AbstractCCCorrection() {
		;
	}

	virtual ~AbstractCCCorrection() {
		;
	}

	virtual double CCEnergy(const compchem::AbstractMatrix<double> &orbitals,
	        const compchem::AbstractMatrix<double> &fock,
	        const compchem::AbstractMatrix<double> &teri,
	        int nelectrons) = 0;
};

}

#endif /* PROJECT_5_CC_HPP_ */
