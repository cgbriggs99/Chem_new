/*
 * ccsdt.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: cgbri
 */

#ifndef PROJECT_6_CCSDT_HPP_
#define PROJECT_6_CCSDT_HPP_

#include "../Project-5/cc.hpp"

namespace compchem {
namespace strategies {

class DefaultCCSDTCorrection: public AbstractCCCorrection {
public:
	DefaultCCSDTCorrection() {
		;
	}

	virtual ~DefaultCCSDTCorrection() {
		;
	}

	double CCEnergy(const compchem::AbstractMatrix<double> &orbitals,
	                const compchem::AbstractMatrix<double> &fock,
	                const compchem::AbstractMatrix<double> &teri,
	                int nelectrons);
};

}
}

#endif /* PROJECT_6_CCSDT_HPP_ */
