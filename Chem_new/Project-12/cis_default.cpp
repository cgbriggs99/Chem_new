/*
 * cis_default.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: connor
 */

#ifndef __CIS_DEFAULT_CPP__
#define __CIS_DEFAULT_CPP__
#include "cis_default.hpp"

template<typename _Eigs>
compchem::Matrix<double> &compchem::strategies::DefaultCISEnergyStrategy<_Eigs>::CIEnergy(
        const compchem::AbstractMatrix<double> &sofock,
        const compchem::AbstractMatrix<double> &sotei, int occupied) {

	compchem::Matrix<double> *ham = new compchem::Matrix<double>(
	        {(sofock.getShape(0) - occupied) * occupied, (sofock.getShape(0)
	                - occupied) * occupied});

	//Calculate the new Hamiltonian.
	for(int i = 0, ia = 0; i < occupied; i++) {
		for(int a = occupied; a < sofock.getShape(0); a++, ia++) {
			for(int j = 0, jb = 0; j < occupied; j++) {
				for(int b = occupied; b < sofock.getShape(0); b++, jb++) {
					ham->setEntry(
					        sofock.getEntry(a, b) * ((i == j)? 1: 0)
					                - sofock.getEntry(i, j) * ((a == b)? 1: 0)
					                + sotei.getEntry(a, j, i, b), ia, jb);
				}
			}
		}
	}

	compchem::Matrix<double> &out = this->eigs->eigenvals(*ham);
	delete ham;
	return (out);
}

#endif
