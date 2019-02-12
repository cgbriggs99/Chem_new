/*
 * scf.hpp
 *
 *  Created on: Feb 4, 2019
 *      Author: Connor
 */

#ifndef PROJECT_3_SCF_HPP_
#define PROJECT_3_SCF_HPP_

#include "../Base/matrix.hpp"
#include "../Molecule/wavefunction.hpp"
#include "../Molecule/molecule.hpp"
#include <vector>
#include <array>

namespace compchem {

class SCFStrategy {
public:
	SCFStrategy() {
		;
	}
	virtual ~SCFStrategy() = default;

	virtual compchem::AbstractMatrix<double> &findHamiltonian(const compchem::AbstractWavefunction &wf) = 0;
	virtual void runSCF(const compchem::AbstractWavefunction &wf, compchem::AbstractMatrix<double> **mo_fock,
			compchem::AbstractMatrix<double> **lcaomo, compchem::AbstractMatrix<double> **density, double *energy) = 0;
	virtual std::vector<double> &findElectronCharge(const compchem::AbstractMolecule &mol, const compchem::AbstractWavefunction &wf,
			const compchem::AbstractMatrix<double> &density) = 0;
	virtual std::array<double, 3> &findDipole(const compchem::AbstractMolecule &mol, const compchem::AbstractMatrix<double> &density,
			const compchem::AbstractWavefunction &wf) = 0;
};

}

#endif /* PROJECT_3_SCF_HPP_ */
