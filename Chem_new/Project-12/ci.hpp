/*
 * ci.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: connor
 */

#ifndef PROJECT_12_CI_HPP_
#define PROJECT_12_CI_HPP_

#include "../Base/base.hpp"

namespace compchem {

class AbstractCIEnergyStrategy {
public:
	AbstractCIEnergyStrategy() {
		;
	}
	virtual ~AbstractCIEnergyStrategy() {
		;
	}

	virtual compchem::Matrix<double> &CIEnergy(
	        const compchem::AbstractMatrix<double> &sofock,
	        const compchem::AbstractMatrix<double> &sotei, int occupied) = 0;
};

}

#endif /* PROJECT_12_CI_HPP_ */
