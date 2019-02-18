/*
 * ccsd_default.hpp
 *
 *  Created on: Feb 11, 2019
 *      Author: connor
 */

#ifndef PROJECT_5_CCSD_DEFAULT_HPP_
#define PROJECT_5_CCSD_DEFAULT_HPP_

#include "ccsd.hpp"

namespace compchem {
namespace strategies {

class DefaultCCSDCorrection : public AbstractCCSDCorrection {
public:
	DefaultCCSDCorrection() {
		;
	}

	virtual ~DefaultCCSDCorrection() {
		;
	}

	double CCSDEnergy(const compchem::AbstractMatrix<double> &orbitals,
	        const compchem::AbstractMatrix<double> &fock,
	        const compchem::AbstractMatrix<double> &teri,
	        const compchem::AbstractMatrix<double> &energies, int nelectrons,
	        compchem::AbstractMatrix<double> **t1_amps,
	        compchem::AbstractMatrix<double> **t2_amps) override;
};

}
}

#endif /* PROJECT_5_CCSD_DEFAULT_HPP_ */
