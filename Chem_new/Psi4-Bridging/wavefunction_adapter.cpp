/*
 * wavefunction_adapter.cpp
 *
 *  Created on: Feb 20, 2019
 *      Author: connor
 */

#include "wavefunction_adapter.hpp"
#include "matrix_adapter.hpp"
#include "tei_adapter.hpp"
#include "../Molecule/sto3g_basis_set.hpp"
#include "psi4/libmints/basisset.h"

#ifndef __WAVEFUNCTION_ADAPTER_CPP__
#define __WAVEFUNCTION_ADAPTER_CPP__

//compchem::AbstractWavefunction &compchem::WavefunctionAdapter::convert(
//        const psi::Wavefunction &wf, const psi::TwoElectronInt &tei) const {
//	compchem::strategies::DefaultWavefunction<compchem::strategies::STO3GBasisSet> *out =
//	        new compchem::strategies::DefaultWavefunction<compchem::strategies::STO3GBasisSet>(wf.nso(),
//	                wf.nalpha() + wf.nbeta());
//	const compchem::Psi4MatrixAdapter *matrix= compchem::Psi4MatrixAdapter::getSingleton();
//	const compchem::TEIAdapter *teiadapt = compchem::TEIAdapter::getSingleton();
//	out->setS(&matrix->convert(wf.S()));
//	out->setTEI(&teiadapt->convert(tei));
//	out->setHam(&matrix->convert(wf.H()));
//
//	return (*out);
//}

compchem::AbstractWavefunction &compchem::WavefunctionAdapter::convert(psi::MintsHelper &helper, double enuc) const {
	compchem::strategies::DefaultWavefunction<compchem::strategies::STO3GBasisSet> *out =
			new compchem::strategies::DefaultWavefunction<compchem::strategies::STO3GBasisSet>(helper.nbf(),
					helper.basisset().get()->n_ecp_core());
	const compchem::Psi4MatrixAdapter *matrix= compchem::Psi4MatrixAdapter::getSingleton();
	const compchem::TEIAdapter *teiadapt = compchem::TEIAdapter::getSingleton();

	std::vector<psi::SharedMatrix> dipole = helper.ao_dipole();

	out->setMuX(&matrix->convert(dipole[0]));
	out->setMuY(&matrix->convert(dipole[1]));
	out->setMuZ(&matrix->convert(dipole[2]));
	out->setS(&matrix->convert(helper.ao_overlap()));
	out->setT(&matrix->convert(helper.ao_potential()));
	out->setV(&matrix->convert(helper.ao_kinetic()));
	out->setTEI(&teiadapt->convert(helper.ao_eri()));
	out->setEnuc(enuc);
	return (*out);
}

#endif
