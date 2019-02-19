/*
 * cis_default.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: connor
 */

#ifndef PROJECT_12_CIS_DEFAULT_HPP_
#define PROJECT_12_CIS_DEFAULT_HPP_

#include "ci.hpp"

namespace compchem {
namespace strategies {

template<typename _Eigs>
class DefaultCISEnergyStrategy : public AbstractCIEnergyStrategy {
protected:
	_Eigs *eigs;
public:
	DefaultCISEnergyStrategy() {
		eigs = new _Eigs();
	}

	virtual ~DefaultCISEnergyStrategy() {
		delete eigs;
	}

	compchem::Matrix<double> &CIEnergy(
	        const compchem::AbstractMatrix<double> &sofock,
	        const compchem::AbstractMatrix<double> &sotei, int occupied);
};

}
}

#include "cis_default.cpp"

#endif /* PROJECT_12_CIS_DEFAULT_HPP_ */
