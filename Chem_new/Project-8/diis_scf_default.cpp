/*
 * diis_scf_default.cpp
 *
 *  Created on: Feb 13, 2019
 *      Author: connor
 */

#include "diis_scf.hpp"
#include <queue>

#ifndef __DIIS_SCF_DEFAULT_CPP__
#define __DIIS_SCF_DEFAULT_CPP__

template<typename _T, typename _U, int _max_depth>
void compchem::strategies::SCF_DIISStrategy<_T, _U, _max_depth>::runSCF(
        const compchem::AbstractWavefunction &wf,
        compchem::AbstractMatrix<double> **mo_fock,
        compchem::AbstractMatrix<double> **lcaomo,
        compchem::AbstractMatrix<double> **density,
        compchem::AbstractMatrix<double> **eigs, double *energy) {

	std::queue<compchem::Matrix<double> *> errmats, fockmats;

	compchem::Matrix<double> *s_half, *s_half_t, *fock = nullptr, *fock_ao,
	        *c_prime = nullptr, *c = nullptr, *dens = nullptr, *last_dens =
	                nullptr, *work1, *work2, *work3, *eout = nullptr;
	compchem::Matrix<double> *hamiltonian =
	        (compchem::Matrix<double> *) &findHamiltonian(wf);
	double etotal = 0, etot_last = 1, rms = 0, rms_last = 1;

	/*
	 * Calculate S^-1/2
	 */
	this->eigs->eigen_all(wf.s(), &s_half, &work1, nullptr);

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

	while(fabs(etotal - etot_last) > 0.000000001
	        || fabs(rms - rms_last) > 0.000001) {
		if(last_dens != nullptr) {
			delete last_dens;
		}

		if(c != nullptr)
			delete c;
		if(c_prime != nullptr)
			delete c_prime;
		if(eout != nullptr) {
			delete eout;
		}

		last_dens = dens;

		etot_last = etotal;
		rms_last = rms;

		//Do DIIS.
		if(dens != nullptr) {
			work1 = (compchem::Matrix<double> *) &matarit->mult(*fock, *dens);
			work2 = (compchem::Matrix<double> *) &matarit->mult(*work1, wf.s());
			delete work1;
			work1 = (compchem::Matrix<double> *) &matarit->mult(wf.s(), *dens);
			work3 = (compchem::Matrix<double> *) &matarit->mult(*work1, *fock);
			delete work1;
			work1 = (compchem::Matrix<double> *) &matarit->sub(*work2, *work3);
			delete work2;
			delete work3;
			errmats.push(work1);
			fockmats.push(fock);

			//If work1 is deleted, the stored matrix will be deleted as well. Nullify to prevent this.
			work1 = nullptr;

			//Keep the queue the right size.
			if(errmats.size() > (_max_depth)) {
				work1 = errmats.front();
				errmats.pop();
				delete work1;
			}
			if(fockmats.size() > (_max_depth)) {
				work1 = fockmats.front();
				errmats.pop();
				delete work1;
			}

			//Build new minimization array.
			work1 = new compchem::Matrix<double>({_max_depth + 1, _max_depth + 1});
			work1->setEntry(0, _max_depth, _max_depth);
			work2 = new compchem::Matrix<double>({_max_depth + 1, 1});

			for(int i = 0; i < _max_depth; i++) {
				for(int j = 0; j < i; j++) {
					double sum = 0;
					compchem::AbstractMatrix<double> *e1, *e2;
					e1 = mats[i];
					e2 = mats[j];
					for(int k = 0; k < wf.getSize(); k++) {
						for(int l = 0; l < wf.getSize(); l++) {
							sum += e1->getEntry(k, l) * e2.getEntry(k, l);
						}
					}
					work1->setEntry(sum, i, j);
					work1->setEntry(sum, j, i);
				}
			}

			for(int i = 0; i < _max_depth; i++) {
				work1->setEntry(-1, i, _max_depth);
				work1->setEntry(-1, _max_depth, i);
			}

			//Solve.
			matarit->solve(*work1, *work2);

			delete work1;

			for(int i = 0; i < _max_depth - 1; i++) {

			}
		}

		//Find mo Fock matrix
		work1 = (compchem::Matrix<double> *) &matarit->mult(*s_half_t,
		        *fock_ao);
		fock = (compchem::Matrix<double> *) &matarit->mult(*work1, *s_half);
		delete work1;

		//Diagonalize Fock matrix
		this->eigs->eigen_all(*fock, &eout, &c_prime, nullptr);

		//Find the lcao.
		c = (compchem::Matrix<double> *) &matarit->mult(*s_half, *c_prime);

		//Calculate the density
		dens = new compchem::Matrix<double>( {wf.getSize(), wf.getSize()});
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
				etotal +=
				        dens->getEntry(i, j)
				                * (hamiltonian->getEntry(i, j)
				                        + fock_ao->getEntry(i, j));
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
						sum +=
						        dens->getEntry(k, l)
						                * (2
						                        * wf.two_electron().getEntry(i,
						                                j, k, l)
						                        - wf.two_electron().getEntry(i,
						                                k, j, l));
					}
				}
				fock_ao->setEntry(sum, i, j);
			}
		}
		continue;
	}
	if(mo_fock != nullptr) {
		*mo_fock = fock;
	} else {
		delete fock;
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
	delete fock_ao;
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

#endif

