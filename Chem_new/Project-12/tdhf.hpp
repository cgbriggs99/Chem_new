/*
 * tdhf.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: connor
 */

#ifndef PROJECT_12_TDHF_HPP_
#define PROJECT_12_TDHF_HPP_
#include "../Base/base.hpp"

namespace compchem {

class TDHFEnergyStrategy {
public:
	TDHFEnergyStrategy() {
		;
	}

	virtual ~TDHFEnergyStrategy() {
		;
	}

	virtual compchem::AbstractMatrix<double> &TDHFEnergies(
	        const compchem::AbstractMatrix<double> &sofock,
	        const compchem::AbstractMatrix<double> &sotei, int occupied) = 0;
};

}

#endif /* PROJECT_12_TDHF_HPP_ */
