/*
 * spin_adapt_default.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: connor
 */

#include "spin_adapt_default.hpp"

compchem::AbstractMatrix<double> &compchem::strategies::DefaultSpinAdaptationStrategy::spinAdaptTEI(
        const compchem::AbstractMatrix<double> &tei,
        const compchem::AbstractMatrix<double> &orbs) {
	compchem::Matrix<double> *out = new compchem::Matrix<double>(
	        {2 * tei.getShape(0), 2 * tei.getShape(0), 2 * tei.getShape(0), 2
	                * tei.getShape(0)}), *motei =
	        new compchem::Matrix<double>(
	                {tei.getShape(0), tei.getShape(0), tei.getShape(0), tei
	                        .getShape(0)});

	for(int p = 0; p < tei.getShape(0); p++) {
		for(int q = 0; q < tei.getShape(0); q++) {
			for(int r = 0; r < tei.getShape(0); r++) {
				for(int s = 0; s < tei.getShape(0); s++) {
					double sum = 0;
					for(int mu = 0; mu < tei.getShape(0); mu++) {
						double sum1 = 0;
						for(int nu = 0; nu < tei.getShape(0); nu++) {
							double sum2 = 0;
							for(int lam = 0; lam < tei.getShape(0); lam++) {
								double sum3 = 0;
								for(int sig = 0; sig < tei.getShape(0); sig++) {
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

	for(int p = 0; p < 2 * tei.getShape(0); p++) {
		for(int q = 0; q < 2 * tei.getShape(0); q++) {
			for(int r = 0; r < 2 * tei.getShape(0); r++) {
				for(int s = 0; s < 2 * tei.getShape(0); s++) {
					out->setEntry(
					        motei->getEntry(p / 2, r / 2, q / 2, s / 2)
					                * ((p % 2 == r % 2)? 1: 0)
					                * ((q % 2 == s % 2)? 1: 0)
					                - motei->getEntry(p / 2, s / 2, q / 2,
					                        r / 2) * ((p % 2 == s % 2)? 1: 0)
					                        * ((q % 2 == r % 2)? 1: 0), p, q, r,
					        s);
				}
			}
		}
	}

	delete motei;
	return (*out);
}

compchem::AbstractMatrix<double> &compchem::strategies::DefaultSpinAdaptationStrategy::spinAdapt1elec(
        const compchem::AbstractMatrix<double> &arr,
        const compchem::AbstractMatrix<double> &orbitals) {
	compchem::Matrix<double> *out = new compchem::Matrix<double>(
	        {2 * arr.getShape(0), 2 * arr.getShape(0)}), *moarr =
	        new compchem::Matrix<double>( {arr.getShape(0), arr.getShape(0)});

	for(int i = 0; i < arr.getShape(0); i++) {
		for(int j = 0; j < arr.getShape(0); j++) {
			double sum1 = 0;
			for(int k = 0; k < arr.getShape(0); k++) {
				double sum2 = 0;
				for(int l = 0; l < arr.getShape(0); l++) {
					sum2 += orbitals.getEntry(l, j) * arr.getEntry(k, l);
				}
				sum1 += orbitals.getEntry(k, i) * sum2;
			}
			moarr->setEntry(sum1, i, j);
		}
	}

	for(int i = 0; i < 2 * arr.getShape(0); i++) {
		for(int j = 0; j < 2 * arr.getShape(0); j++) {
			out->setEntry(
			        moarr->getEntry(i / 2, j / 2) * ((i % 2 == j % 2)? 1: 0),
			        i, j);
		}
	}
	delete moarr;

	return (*out);
}
