/*
 * ccsd.hpp
 *
 *  Created on: Feb 11, 2019
 *      Author: connor
 */

#ifndef PROJECT_5_CCSD_HPP_
#define PROJECT_5_CCSD_HPP_

#include "../Base/base.hpp"

namespace compchem {

class AbstractCCSDCorrection {
public:
	AbstractCCSDCorrection() {
		;
	}

	virtual ~AbstractCCSDCorrection() {
		;
	}

	virtual double CCSDEnergy(const compchem::AbstractMatrix<double> &orbitals,
			const compchem::AbstractMatrix<double> &fock,
			const compchem::AbstractMatrix<double> &teri,
			const compchem::AbstractMatrix<double> &energies, const compchem::AbstractMatrix<double> &hamiltonian,
			int nelectrons) = 0;
};

}

#endif /* PROJECT_5_CCSD_HPP_ */
