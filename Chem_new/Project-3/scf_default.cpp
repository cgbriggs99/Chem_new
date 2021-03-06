/*
 * scf_default.cpp
 *
 *  Created on: Feb 4, 2019
 *      Author: Connor
 */

#ifndef __SCF_DEFAULT_CPP__
#define __SCF_DEFAULT_CPP__

#include "scf_default.hpp"

template<typename _E, typename _M>
compchem::AbstractMatrix<double> &compchem::strategies::DefaultSCFStrategy<_E, _M>::findHamiltonian(const compchem::AbstractWavefunction &wf) {
	compchem::Matrix<double> *out = new compchem::Matrix<double>({wf.getSize(), wf.getSize()});
	const compchem::AbstractMatrix<double> *t = &wf.t(), *v = &wf.v();

	for(int i = 0; i < wf.getSize(); i++) {
		for(int j = 0; j < wf.getSize(); j++) {
			out->setEntry(t->getEntry(i, j) + v->getEntry(i, j), i, j);
		}
	}
	return (*out);
}

template<typename _E, typename _M>
void compchem::strategies::DefaultSCFStrategy<_E, _M>::runSCF(const compchem::AbstractWavefunction &wf, compchem::AbstractMatrix<double> **mo_fock,
		compchem::AbstractMatrix<double> **lcaomo, compchem::AbstractMatrix<double> **density,
		compchem::AbstractMatrix<double> **_eigs, double *energy) {

	compchem::Matrix<double> *s_half, *s_half_t, *fock = nullptr, *fock_ao, *c_prime = nullptr, *c = nullptr, *dens = nullptr,
			*last_dens = nullptr, *work1, *work2, *eout = nullptr;
	compchem::Matrix<double> *hamiltonian = (compchem::Matrix<double> *) &findHamiltonian(wf);
	double etotal = 0, etot_last = 1, rms = 0, rms_last = 1;

	/*
	 * Calculate S^-1/2
	 */
	eigs->eigen_all(wf.s(), &s_half, &work1, nullptr);

	for(int i = 0; i < s_half->getShape(0); i++) {
		s_half->setEntry(1.0 / sqrt(s_half->getEntry(i)), i);
	}

	work2 = (compchem::Matrix<double> *) &matarit->inverse(*work1);

	for(int i = 0; i < work1->getShape(0); i++) {
		for(int j = 0; j < work1->getShape(1); j++) {
			work1->setEntry(s_half->getEntry(j) * work1->getEntry(i, j), i, j);
		}
	}

	delete s_half;
	s_half = (compchem::Matrix<double> *) &matarit->mult(*work1, *work2);
	s_half_t = (compchem::Matrix<double> *) &matarit->transpose(*s_half);
	delete work2;
	delete work1;

	//Set up initial Fock
	fock_ao = (compchem::Matrix<double> *) &findHamiltonian(wf);

	while(fabs(etotal - etot_last) > 0.000000001 || fabs(rms - rms_last) > 0.000001) {
		if(last_dens != nullptr) {
			delete last_dens;
		}

		if(c != nullptr)
			delete c;
		if(c_prime != nullptr)
			delete c_prime;
		if(fock != nullptr)
			delete fock;
		if(eout != nullptr) {
			delete eout;
		}

		last_dens = dens;

		etot_last = etotal;
		rms_last = rms;

		//Find mo Fock matrix
		work1 = (compchem::Matrix<double> *) &matarit->mult(*s_half_t, *fock_ao);
		fock = (compchem::Matrix<double> *) &matarit->mult(*work1, *s_half);
		delete work1;


		//Diagonalize Fock matrix
		eigs->eigen_all(*fock, &eout, &c_prime, nullptr);

		//Find the lcao.
		c = (compchem::Matrix<double> *) &matarit->mult(*s_half, *c_prime);

		//Calculate the density
		dens = new compchem::Matrix<double>({wf.getSize(), wf.getSize()});
		for(int i = 0; i < wf.getSize(); i++) {
			for(int j = 0; j < wf.getSize(); j++) {
				double sum = 0;
				for(int k = 0; k < wf.nelectron() / 2; k++) {
					sum += c->getEntry(i, k) * c->getEntry(j, k);
				}
				dens->setEntry(sum, i, j);
			}
		}

		//Calculate the energy.
		etotal = wf.enuc();
		for(int i = 0; i < wf.getSize(); i++) {
			for(int j = 0; j < wf.getSize(); j++) {
				etotal += dens->getEntry(i, j) * (hamiltonian->getEntry(i, j) + fock_ao->getEntry(i, j));
			}
		}

		if(last_dens != nullptr) {
			rms = 0;
			for(int i = 0; i < wf.getSize(); i++) {
				for(int j = 0; j < wf.getSize(); j++) {
					rms += dens->getEntry(i, j) - last_dens->getEntry(i, j);
				}
			}
		} else {
			rms = 0;
			rms_last = 1;
		}

		for(int i = 0; i < wf.getSize(); i++) {
			for(int j = 0; j < wf.getSize(); j++) {
				double sum = hamiltonian->getEntry(i, j);
				for(int k = 0; k < wf.getSize(); k++) {
					for(int l = 0; l < wf.getSize(); l++) {
						sum += dens->getEntry(k, l) * (2 * wf.two_electron().getEntry(i, j, k, l) - wf.two_electron().getEntry(i, k, j, l));
					}
				}
				if(sum > 100) {
					;
				}
				fock_ao->setEntry(sum, i, j);
			}
		}
		continue;
	}
	if(mo_fock != nullptr) {
		*mo_fock = fock_ao;
	} else {
		delete fock_ao;
	}
	if(lcaomo != nullptr) {
		*lcaomo = c;
	} else {
		delete c;
	}
	if(density != nullptr) {
		*density = dens;
	} else {
		delete dens;
	}
	if(_eigs != nullptr) {
		*_eigs = eout;
	} else {
		delete eout;
	}
	delete s_half;
	delete s_half_t;
	delete fock;
	delete hamiltonian;
	if(c_prime != nullptr)
		delete c_prime;
	if(last_dens != nullptr) {
		delete last_dens;
	}
	if(energy != nullptr) {
		*energy = etotal;
	}
}

template<typename _E, typename _M>
std::vector<double> &compchem::strategies::DefaultSCFStrategy<_E, _M>::findElectronCharge(const compchem::AbstractMolecule &mol,
		const compchem::AbstractWavefunction &wf,
		const compchem::AbstractMatrix<double> &density) {
	compchem::AbstractMatrix<double> *ds = &this->matarit->mult(density, wf.s());
	std::vector<double> *out = new std::vector<double>(mol.natom());
	int ind = 0;

	for(int i = 0; i < mol.natom(); i++) {
		double sum = mol.fZ(i);
		for(int j = ind; j < ind + wf.norbitals((int) mol.fZ(i)); j++) {
			sum -= 2 * (ds->getEntry(j, j));
		}
		ind += wf.norbitals((int) mol.fZ(i));
		out->at(i) = sum;
	}
	delete ds;
	return (*out);
}

template<typename _E, typename _M>
std::array<double, 3> &compchem::strategies::DefaultSCFStrategy<_E, _M>::findDipole(const compchem::AbstractMolecule &mol, const compchem::AbstractMatrix<double> &density,
			const compchem::AbstractWavefunction &wf) {
	std::array<double, 3> *out = new std::array<double, 3>();
	out->at(0) = 0;
	out->at(1) = 0;
	out->at(2) = 0;
	for(int i = 0; i < wf.getSize(); i++) {
		for(int j = 0; j < wf.getSize(); j++) {
			out->at(0) += 2 * density.getEntry(i, j) * wf.mux().getEntry(i, j);
			out->at(1) += 2 * density.getEntry(i, j) * wf.muy().getEntry(i, j);
			out->at(2) += 2 * density.getEntry(i, j) * wf.muz().getEntry(i, j);
		}
	}

	for(int i = 0; i < mol.natom(); i++) {
		out->at(0) += mol.fZ(i) * mol.x(i);
		out->at(1) += mol.fZ(i) * mol.y(i);
		out->at(2) += mol.fZ(i) * mol.z(i);
	}
	return (*out);
}

#endif
