/*
 * ccsd_default.hpp
 *
 *  Created on: Feb 11, 2019
 *      Author: connor
 */

#ifndef PROJECT_5_CCSD_DEFAULT_HPP_
#define PROJECT_5_CCSD_DEFAULT_HPP_

#include "cc.hpp"

namespace compchem {
namespace strategies {

class DefaultCCSDCorrection : public AbstractCCCorrection {
public:
	DefaultCCSDCorrection() {
		;
	}

	virtual ~DefaultCCSDCorrection() {
		;
	}

	double CCEnergy(const compchem::AbstractMatrix<double> &orbitals,
	        const compchem::AbstractMatrix<double> &fock,
	        const compchem::AbstractMatrix<double> &teri, int nelectrons) override;
};

}
}

#endif /* PROJECT_5_CCSD_DEFAULT_HPP_ */
