/*
 * mp2_default.cpp
 *
 *  Created on: Feb 11, 2019
 *      Author: connor
 */

#ifndef __MP2_DEFAULT_CPP__
#define __MP2_DEFAULT_CPP__

#include "mp2_default.hpp"
#include "../Base/base.hpp"

template<typename _T>
double compchem::strategies::DefaultMP2Strategy<_T>::mp2Energy(
		const compchem::AbstractMolecule &mol,
		const compchem::AbstractMatrix<double> &tei,
		const compchem::AbstractMatrix<double> &orbs,
		const compchem::AbstractMatrix<double> &eigs) {
	int orbitals = basis_set->norbitals(mol);
	compchem::Matrix<double> *motei = new compchem::Matrix<double>( { orbitals,
			orbitals, orbitals, orbitals });

	for (int p = 0; p < orbitals; p++) {
		for (int q = 0; q < orbitals; q++) {
			for (int r = 0; r < orbitals; r++) {
				for (int s = 0; s < orbitals; s++) {
					double sum = 0;
					for (int mu = 0; mu < orbitals; mu++) {
						double sum1 = 0;
						for (int nu = 0; nu < orbitals; nu++) {
							double sum2 = 0;
							for (int lam = 0; lam < orbitals; lam++) {
								double sum3 = 0;
								for (int sig = 0; sig < orbitals; sig++) {
									sum3 += orbs.getEntry(sig, s)
											* tei.getEntry(mu, nu, lam, sig);
								}
								sum2 += orbs.getEntry(lam, r) * sum3;
							}
							sum1 += orbs.getEntry(nu, q) * sum2;
						}
						sum += orbs.getEntry(mu, p) * sum1;
					}
					motei->setEntry(sum, p, q, r, s);
				}
			}
		}
	}

	double energy = 0;
	for (int i = 0; i < mol.nelectron() / 2; i++) {
		for (int j = 0; j < mol.nelectron() / 2; j++) {
			for (int a = mol.nelectron() / 2; a < orbitals; a++) {
				for (int b = mol.nelectron() / 2; b < orbitals; b++) {
					energy += (motei->getEntry(i, a, j, b)
							* (2 * motei->getEntry(i, a, j, b)
									- motei->getEntry(i, b, j, a)))
							/ (eigs.getEntry(i) + eigs.getEntry(j)
									- eigs.getEntry(a) - eigs.getEntry(b));
				}
			}
		}
	}
	delete motei;
	return (energy);
}

#endif
