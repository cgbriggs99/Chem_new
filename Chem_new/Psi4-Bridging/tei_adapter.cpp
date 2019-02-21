/*
 * tei_adapter.cpp
 *
 *  Created on: Feb 20, 2019
 *      Author: connor
 */

#ifndef __TEI_ADAPTER_CPP__
#define __TEI_ADAPTER_CPP__

#include "tei_adapter.hpp"
#include <math.h>

compchem::strategies::TEIMatrix<double> &compchem::TEIAdapter::convert(
        psi::SharedMatrix tei) const {
	compchem::strategies::TEIMatrix<double> *out =
	        new compchem::strategies::TEIMatrix<double>(
	                (int) sqrt((double) tei.get()->nrow()));

	for(int i = 0; i < out->getShape(0); i++) {
		for(int j = 0; j <= i; j++) {
			for(int k = 0; k <= i; k++) {
				for(int l = 0; l <= (i == k)? j: k; l++) {
					out->setEntry(tei.get()->get(i * out->getShape(0) + j, k * out->getShape(0) + l), i, j, k, l);
				}
			}
		}
	}

	return (*out);
}

#endif
