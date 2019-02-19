/*
 * tdhf_rpa2.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: connor
 */

#ifndef __TDHF_RPA2_CPP__
#define __TDHF_RPA2_CPP__

#include "tdhf_rpa2.hpp"
#include <math.h>

template<typename _Eigs, typename _Matarit>
compchem::AbstractMatrix<double> &compchem::strategies::EERPATDHFEnergyStrategy<
        _Eigs, _Matarit>::TDHFEnergies(const compchem::AbstractMatrix<double> &sofock,
        const compchem::AbstractMatrix<double> &sotei, int occupied) {
	int offset = occupied * (sofock.getShape(0) - occupied);
	compchem::Matrix<double> *mat_a = new compchem::Matrix<double>( {offset,
	        offset}), *mat_b = new compchem::Matrix<double>( {offset, offset});

	for(int i = 0, ia = 0; i < occupied; i++) {
		for(int a = occupied; a < sofock.getShape(0); a++, ia++) {
			for(int j = 0, jb = 0; j < occupied; j++) {
				for(int b = occupied; b < sofock.getShape(0); b++, jb++) {
					mat_a->setEntry(sofock.getEntry(a, b) * ((i == j)? 1: 0)
					                - sofock.getEntry(i, j) * ((a == b)? 1: 0)
					                + sotei.getEntry(a, j, i, b), ia, jb);
					mat_b->setEntry(sotei.getEntry(a, b, i, j), ia, jb);
				}
			}
		}
	}

	compchem::AbstractMatrix<double> *apb = &mat->add(*mat_a, *mat_b), *amb = &mat->subtract(*mat_a, *mat_b), *prod = &mat->mult(*apb, *amb);
	compchem::Matrix<double> &out = eigs->eigenvals(*prod);

	for(int i = 0; i < out.getShape(0); i++) {
		out.setEntry(sqrt(out.getEntry(i)), i);
	}

	delete mat_a;
	delete mat_b;
	delete apb;
	delete amb;
	delete prod;

	return (out);
}

#endif
