/*
 * spin_adapt.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: connor
 */

#ifndef SPIN_ADAPT_HPP_
#define SPIN_ADAPT_HPP_

#include "../Base/base.hpp"

namespace compchem {

class AbstractSpinAdaptationStrategy {
public:
	AbstractSpinAdaptationStrategy() {
		;
	}

	virtual ~AbstractSpinAdaptationStrategy() {
		;
	}

	virtual compchem::AbstractMatrix<double> &spinAdaptTEI(
	        const compchem::AbstractMatrix<double> &tei,
	        const compchem::AbstractMatrix<double> &orbitals) = 0;

	virtual compchem::AbstractMatrix<double> &spinAdapt1elec(
	        const compchem::AbstractMatrix<double> &arr,
	        const compchem::AbstractMatrix<double> &orbitals) = 0;
};

}

#endif /* SPIN_ADAPT_HPP_ */
