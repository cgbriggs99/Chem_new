/*
 * tdhf_rpa1.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: connor
 */

#ifndef PROJECT_12_TDHF_RPA1_HPP_
#define PROJECT_12_TDHF_RPA1_HPP_

#include "tdhf.hpp"

namespace compchem {
namespace strategies {

template<typename _Eigs>
class RPATDHFEnergyStrategy : public TDHFEnergyStrategy {
protected:
	_Eigs *eigs;
public:
	RPATDHFEnergyStrategy() {
		eigs = new _Eigs();
	}

	virtual ~RPATDHFEnergyStrategy() {
		delete eigs;
	}

	compchem::AbstractMatrix<double> &TDHFEnergies(
		        const compchem::AbstractMatrix<double> &sofock,
		        const compchem::AbstractMatrix<double> &sotei, int occupied) override;

};


}
}

#include "tdhf_rpa1.cpp"

#endif /* PROJECT_12_TDHF_RPA1_HPP_ */
