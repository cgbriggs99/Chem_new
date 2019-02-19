/*
 * tdhf_rpa2.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: connor
 */

#ifndef PROJECT_12_TDHF_RPA2_HPP_
#define PROJECT_12_TDHF_RPA2_HPP_

#include "tdhf.hpp"

namespace compchem {
namespace strategies {

template<typename _Eigs, typename _Matarit>
class EERPATDHFEnergyStrategy : public TDHFEnergyStrategy {
private:
	_Eigs *eigs;
	_Matarit *mat;
public:
	EERPATDHFEnergyStrategy() {
		eigs = new _Eigs();
		mat = new _Matarit();
	}

	virtual ~EERPATDHFEnergyStrategy() {
		delete eigs;
		delete mat;
	}

	compchem::AbstractMatrix<double> &TDHFEnergies(
	        const compchem::AbstractMatrix<double> &sofock,
	        const compchem::AbstractMatrix<double> &sotei, int occupied)
	                override;

};

}
}

#include "tdhf_rpa2.cpp"

#endif /* PROJECT_12_TDHF_RPA2_HPP_ */
