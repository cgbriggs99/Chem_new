/*
 * scf_default.hpp
 *
 *  Created on: Feb 4, 2019
 *      Author: Connor
 */

#ifndef PROJECT_3_SCF_DEFAULT_HPP_
#define PROJECT_3_SCF_DEFAULT_HPP_

#include "scf.hpp"
#include "../Base/eigenvalues_default.hpp"
#include "../Base/matrix_arithmetic_default.hpp"

namespace compchem {
namespace strategies {

template<typename _Eigen, typename _Matarit>
class DefaultSCFStrategy : public SCFStrategy {
protected:
	_Eigen *eigs;
	_Matarit *matarit;
public:
	DefaultSCFStrategy() {
		eigs = new _Eigen();
		matarit = new _Matarit();
	}

	virtual ~DefaultSCFStrategy() {
		delete eigs;
		delete matarit;
	}

	compchem::AbstractMatrix<double> &findHamiltonian(const compchem::AbstractWavefunction &wf);

	void runSCF(const compchem::AbstractWavefunction &wf, compchem::AbstractMatrix<double> **mo_fock,
			compchem::AbstractMatrix<double> **lcaomo, compchem::AbstractMatrix<double> **density,
			compchem::AbstractMatrix<double> **eigs, double *energy);

	std::vector<double> &findElectronCharge(const compchem::AbstractMolecule &mol, const compchem::AbstractWavefunction &wf,
			const compchem::AbstractMatrix<double> &density);

	std::array<double, 3> &findDipole(const compchem::AbstractMolecule &mol, const compchem::AbstractMatrix<double> &density,
			const compchem::AbstractWavefunction &wf);
};

}
}

#include "scf_default.cpp"

#endif /* PROJECT_3_SCF_DEFAULT_HPP_ */
