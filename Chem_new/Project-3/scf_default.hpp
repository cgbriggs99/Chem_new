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

class DefaultSCFStrategy : public SCFStrategy {
private:
	compchem::EigenvalueStrategy<double> *eigs;
	compchem::MatrixArithmeticStrategy<double> *matarit;
public:
	DefaultSCFStrategy() {
		eigs = new compchem::strategies::LapackEigenvalues<double>();
		matarit = new compchem::strategies::DefaultMatrixArithmeticStrategy<double>();
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



#endif /* PROJECT_3_SCF_DEFAULT_HPP_ */
