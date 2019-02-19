/*
 * spin_adapt_default.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: connor
 */

#ifndef SPIN_ADAPT_DEFAULT_HPP_
#define SPIN_ADAPT_DEFAULT_HPP_

#include "spin_adapt.hpp"

namespace compchem {

namespace strategies {

class DefaultSpinAdaptationStrategy : public AbstractSpinAdaptationStrategy {
public:
	DefaultSpinAdaptationStrategy() {
		;
	}

	virtual ~DefaultSpinAdaptationStrategy() {
		;
	}

	compchem::AbstractMatrix<double> &spinAdaptTEI(
			const compchem::AbstractMatrix<double> &tei,
			const compchem::AbstractMatrix<double> &orbitals);

	compchem::AbstractMatrix<double> &spinAdapt1elec(
			const compchem::AbstractMatrix<double> &arr,
			const compchem::AbstractMatrix<double> &orbitals);


};

}

}



#endif /* SPIN_ADAPT_DEFAULT_HPP_ */
