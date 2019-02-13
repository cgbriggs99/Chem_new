/*
 * diis_scf.hpp
 *
 *  Created on: Feb 13, 2019
 *      Author: connor
 */

#ifndef PROJECT_8_DIIS_SCF_HPP_
#define PROJECT_8_DIIS_SCF_HPP_

#include "../Project-3/scf_default.hpp"

namespace compchem {
namespace strategies {

template<typename _Eigen, typename _Matarit, int max_depth = 6>
class SCF_DIISStrategy : public DefaultSCFStrategy<_Eigen, _Matarit> {
public:
	SCF_DIISStrategy() {
		;
	}

	virtual ~SCF_DIISStrategy() {
		;
	}

	void runSCF(const compchem::AbstractWavefunction &wf, compchem::AbstractMatrix<double> **mo_fock,
			compchem::AbstractMatrix<double> **lcaomo, compchem::AbstractMatrix<double> **density,
			compchem::AbstractMatrix<double> **eigs, double *energy) override;
};

}
}



#endif /* PROJECT_8_DIIS_SCF_HPP_ */
