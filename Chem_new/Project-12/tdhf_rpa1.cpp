/*
 * tdhf_rpa1.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: connor
 */

#ifndef __TDHF_RPA1_CPP__
#define __TDHF_RPA1_CPP__

#include "tdhf_rpa1.hpp"

template<typename _Eigs>
compchem::AbstractMatrix<double> &compchem::strategies::RPATDHFEnergyStrategy<
        _Eigs>::TDHFEnergies(const compchem::AbstractMatrix<double> &sofock,
        const compchem::AbstractMatrix<double> &sotei, int occupied) {
	compchem::Matrix<double> *mat = new compchem::Matrix<double>(
	        {2 * occupied * (sofock.getShape(0) - occupied), 2 * occupied
	                * (sofock.getShape(0) - occupied)});
	int offset = occupied * (sofock.getShape(0) - occupied);
	for(int i = 0, ia = 0; i < occupied; i++) {
		for(int a = occupied; a < sofock.getShape(0); a++, ia++) {
			for(int j = 0, jb = 0; j < occupied; j++) {
				for(int b = occupied; b < sofock.getShape(0); b++, jb++) {
					mat->setEntry(
					        sofock.getEntry(a, b) * ((i == j)? 1: 0)
					                - sofock.getEntry(i, j) * ((a == b)? 1: 0)
					                + sotei.getEntry(a, j, i, b), ia, jb);
					mat->setEntry(-mat->getEntry(ia, jb), offset + ia, offset + jb);
					mat->setEntry(sotei.getEntry(a, b, i, j), offset + ia, jb);
					mat->setEntry(-sotei.getEntry(a, b, i, j), ia, offset + jb);
				}
			}
		}
	}
	compchem::Matrix<double> &out = eigs->eigenvals(*mat);
	delete mat;
	return (out);
}

#endif
