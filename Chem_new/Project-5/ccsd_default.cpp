/*
 * ccsd_default.cpp
 *
 *  Created on: Feb 11, 2019
 *      Author: connor
 */

#include "ccsd_default.hpp"
#include <cmath>

static compchem::Matrix<double> &calculateSOTEI(
        const compchem::AbstractMatrix<double> &tei,
        const compchem::AbstractMatrix<double> &orbs, int occupied) {
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

static compchem::Matrix<double> &calculateSOFock(
        const compchem::AbstractMatrix<double> &fock,
        const compchem::AbstractMatrix<double> &orbs) {
	compchem::Matrix<double> *out = new compchem::Matrix<double>(
	        {2 * fock.getShape(0), 2 * fock.getShape(0)}), *mofock =
	        new compchem::Matrix<double>( {fock.getShape(0), fock.getShape(0)});

	for(int i = 0; i < fock.getShape(0); i++) {
		for(int j = 0; j < fock.getShape(0); j++) {
			double sum1 = 0;
			for(int k = 0; k < fock.getShape(0); k++) {
				double sum2 = 0;
				for(int l = 0; l < fock.getShape(0); l++) {
					sum2 += orbs.getEntry(l, j) * fock.getEntry(k, l);
				}
				sum1 += orbs.getEntry(k, i) * sum2;
			}
			if(fabs(sum1) < 0.000000000000001) {
				mofock->setEntry(0, i, j);
			} else {
				mofock->setEntry(sum1, i, j);
			}
		}
	}

	for(int i = 0; i < 2 * fock.getShape(0); i++) {
		for(int j = 0; j < 2 * fock.getShape(0); j++) {
			out->setEntry(
			        mofock->getEntry(i / 2, j / 2) * ((i % 2 == j % 2)? 1: 0),
			        i, j);
		}
	}
	delete mofock;

	return (*out);
}

static void calculateFIntermediate(
        const compchem::AbstractMatrix<double> &sofock,
        const compchem::AbstractMatrix<double> &sotei,
        const compchem::AbstractMatrix<double> &t1,
        const compchem::AbstractMatrix<double> &t2, int occupied,
        compchem::AbstractMatrix<double> &out) {
	//Eq. 3: unoccupied-unoccupied.
	for(int a = occupied; a < sofock.getShape(0); a++) {
		for(int e = occupied; e < sofock.getShape(0); e++) {
			double tot = (1 - ((a == e)? 1: 0)) * sofock.getEntry(a, e);
			double sum1 = 0, sum2 = 0, sum3 = 0;

			for(int m = 0; m < occupied; m++) {
				sum1 += sofock.getEntry(m, e) * t1.getEntry(m, a) / 2;
			}

			for(int m = 0; m < occupied; m++) {
				for(int f = occupied; f < sofock.getShape(0); f++) {
					sum2 += t1.getEntry(m, f) * sotei.getEntry(m, a, f, e);
				}
			}

			for(int m = 0; m < occupied; m++) {
				for(int n = 0; n < occupied; n++) {
					for(int f = occupied; f < sofock.getShape(0); f++) {
						sum3 += (t2.getEntry(m, n, a, f)
						        + (t1.getEntry(m, a) * t1.getEntry(n, f)
						                - t1.getEntry(m, f) * t1.getEntry(n, a))
						                / 2) * sotei.getEntry(m, n, e, f) / 2;
					}
				}
			}
			out.setEntry(tot - sum1 + sum2 - sum3, a, e);
		}
	}

	//Eq. 4: occupied-occupied
	for(int m = 0; m < occupied; m++) {
		for(int i = 0; i < occupied; i++) {
			double tot = (1 - ((m == i)? 1: 0)) * sofock.getEntry(m, i), sum1 =
			        0, sum2 = 0, sum3 = 0;

			for(int e = occupied; e < sofock.getShape(0); e++) {
				sum1 += t1.getEntry(i, e) * sofock.getEntry(m, e) / 2;
			}

			for(int e = occupied; e < sofock.getShape(0); e++) {
				for(int n = 0; n < occupied; n++) {
					sum2 += t1.getEntry(n, e) * sotei.getEntry(m, n, i, e);
				}
			}

			for(int e = occupied; e < sofock.getShape(0); e++) {
				for(int n = 0; n < occupied; n++) {
					for(int f = occupied; f < sofock.getShape(0); f++) {
						sum3 += (t2.getEntry(i, n, e, f)
						        + (t1.getEntry(i, e) * t1.getEntry(n, f)
						                - t1.getEntry(i, f) * t1.getEntry(n, e))
						                / 2) * sotei.getEntry(m, n, e, f) / 2;
					}
				}
			}
			out.setEntry(tot + sum1 + sum2 + sum3, m, i);
		}

		//Eq. 5: occupied-unoccupied
		for(int m = 0; m < occupied; m++) {
			for(int e = occupied; e < sofock.getShape(0); e++) {
				double sum = 0;
				for(int n = 0; n < occupied; n++) {
					for(int f = occupied; f < sofock.getShape(0); f++) {
						sum += t1.getEntry(n, f) * sotei.getEntry(m, n, e, f);
					}
				}
				out.setEntry(sofock.getEntry(m, e) + sum, m, e);
			}
		}
	}
}

static void calculateWIntermediate(
        const compchem::AbstractMatrix<double> &sofock,
        const compchem::AbstractMatrix<double> &sotei,
        const compchem::AbstractMatrix<double> &t1,
        const compchem::AbstractMatrix<double> &t2, int occupied,
        compchem::AbstractMatrix<double> &out) {
	//Eq. 6: occupied-occupied
	for(int m = 0; m < occupied; m++) {
		for(int n = 0; n < occupied; n++) {
			for(int i = 0; i < occupied; i++) {
				for(int j = 0; j < occupied; j++) {
					double sum1 = 0, sum2 = 0, sum3 = 0;
					for(int e = occupied; e < sofock.getShape(0); e++) {
						sum1 += t1.getEntry(j, e) * sotei.getEntry(m, n, i, e);
						sum2 += t1.getEntry(i, e) * sotei.getEntry(m, n, j, e);
					}

					for(int e = occupied; e < sofock.getShape(0); e++) {
						for(int f = occupied; f < sofock.getShape(0); f++) {
							sum3 += (t2.getEntry(i, j, e, f)
							        + t1.getEntry(i, e) * t1.getEntry(j, f)
							        - t1.getEntry(i, f) * t1.getEntry(j, e))
							        * sotei.getEntry(m, n, e, f) / 4;
						}
					}
					out.setEntry(
					        sotei.getEntry(m, n, i, j) + sum1 - sum2 + sum3, m,
					        n, i, j);
				}
			}
		}
	}

	//Eq. 7: unoccupied-unoccupied
	for(int a = occupied; a < sofock.getShape(0); a++) {
		for(int b = occupied; b < sofock.getShape(0); b++) {
			for(int e = occupied; e < sofock.getShape(0); e++) {
				for(int f = occupied; f < sofock.getShape(0); f++) {
					double sum1 = 0, sum2 = 0, sum3 = 0;
					for(int m = 0; m < occupied; m++) {
						sum1 += t1.getEntry(m, b) * sotei.getEntry(a, m, e, f);
						sum2 += t1.getEntry(m, a) * sotei.getEntry(b, m, e, f);
					}

					for(int m = 0; m < occupied; m++) {
						for(int n = 0; n < occupied; n++) {
							sum3 += (t2.getEntry(m, n, a, b)
							        + t1.getEntry(m, a) * t1.getEntry(n, b)
							        - t1.getEntry(m, b) * t1.getEntry(n, a))
							        * sotei.getEntry(m, n, e, f) / 4;
						}
					}
					out.setEntry(
					        sotei.getEntry(a, b, e, f) - sum1 + sum2 + sum3, a,
					        b, e, f);
				}
			}
		}
	}

	//Eq. 8: mixed
	for(int m = 0; m < occupied; m++) {
		for(int b = occupied; b < sofock.getShape(0); b++) {
			for(int e = occupied; e < sofock.getShape(0); e++) {
				for(int j = 0; j < occupied; j++) {
					double sum1 = 0, sum2 = 0, sum3 = 0;

					for(int f = occupied; f < sofock.getShape(0); f++) {
						sum1 += t1.getEntry(j, f) * sotei.getEntry(m, b, e, f);
					}

					for(int n = 0; n < sofock.getShape(0); n++) {
						sum2 += t1.getEntry(n, b) * sotei.getEntry(m, n, e, j);
					}

					for(int f = occupied; f < sofock.getShape(0); f++) {
						for(int n = 0; n < occupied; n++) {
							sum3 += (t2.getEntry(j, n, f, b) / 2
							        + t1.getEntry(j, f) * t1.getEntry(n, b))
							        * sotei.getEntry(m, n, e, f);
						}
					}
					out.setEntry(
					        sotei.getEntry(m, b, e, j) + sum1 - sum2 - sum3, m,
					        b, e, j);
				}
			}
		}
	}
}

static void calculateT1(const compchem::AbstractMatrix<double> &sofock,
        const compchem::AbstractMatrix<double> &sotei,
        const compchem::AbstractMatrix<double> &t1,
        const compchem::AbstractMatrix<double> &t2,
        const compchem::AbstractMatrix<double> &f,
        const compchem::AbstractMatrix<double> &w, int occupied,
        compchem::AbstractMatrix<double> &out) {

	//Eq. 1
	for(int i = 0; i < occupied; i++) {
		for(int a = occupied; a < sofock.getShape(0); a++) {
			double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0, sum6 = 0;

			for(int e = occupied; e < sofock.getShape(0); e++) {
				sum1 += t1.getEntry(i, e) * f.getEntry(a, e);
			}

			for(int m = 0; m < occupied; m++) {
				sum2 += t1.getEntry(m, a) * f.getEntry(m, i);
			}

			for(int m = 0; m < occupied; m++) {
				for(int e = occupied; e < sofock.getShape(0); e++) {
					sum3 += t2.getEntry(i, m, a, e) * f.getEntry(m, e);
				}
			}

			for(int n = 0; n < occupied; n++) {
				for(int f = occupied; f < sofock.getShape(0); f++) {
					sum4 += t1.getEntry(n, f) * sotei.getEntry(n, a, i, f);
				}
			}

			for(int m = 0; m < occupied; m++) {
				for(int e = occupied; e < sofock.getShape(0); e++) {
					for(int f = occupied; f < sofock.getShape(0); f++) {
						sum5 += t2.getEntry(i, m, e, f)
						        * sotei.getEntry(m, a, e, f) / 2;
					}
				}
			}

			for(int m = 0; m < occupied; m++) {
				for(int e = occupied; e < sofock.getShape(0); e++) {
					for(int n = 0; n < occupied; n++) {
						sum6 += t2.getEntry(m, n, a, e)
						        * sotei.getEntry(n, m, e, i) / 2;
					}
				}
			}

			if(sofock.getEntry(i, i) - sofock.getEntry(a, a) == 0) {
				throw(new std::exception());
			}
			if(fabs(
			        (sofock.getEntry(i, a) + sum1 - sum2 + sum3 - sum4 - sum5
			                - sum6)
			                / (sofock.getEntry(i, i) - sofock.getEntry(a, a)))
			        >= 10) {
				throw(new std::exception());

			}

			out.setEntry(
			        (sofock.getEntry(i, a) + sum1 - sum2 + sum3 - sum4 - sum5
			                - sum6)
			                / (sofock.getEntry(i, i) - sofock.getEntry(a, a)),
			        i, a);
		}
	}
}

static void calculateT2(const compchem::AbstractMatrix<double> &sofock,
        const compchem::AbstractMatrix<double> &sotei,
        const compchem::AbstractMatrix<double> &t1,
        const compchem::AbstractMatrix<double> &t2,
        const compchem::AbstractMatrix<double> &f,
        const compchem::AbstractMatrix<double> &w, int occupied,
        compchem::AbstractMatrix<double> &out) {
	for(int i = 0; i < occupied; i++) {
		for(int j = 0; j < occupied; j++) {
			for(int a = occupied; a < sofock.getShape(0); a++) {
				for(int b = occupied; b < sofock.getShape(0); b++) {
					double sum1 = 0, sum2 = 0, sum3 = 0, sum4 = 0, sum5 = 0,
					        sum6 = 0, sum7 = 0, sum8 = 0, sum9 = 0, sum10 = 0,
					        sum11 = 0, sum12 = 0, sum13 = 0, sum14 = 0;
					for(int e = occupied; e < sofock.getShape(0); e++) {
						double suma = 0, sumb = 0;
						for(int m = 0; m < occupied; m++) {
							suma += t1.getEntry(m, b) * f.getEntry(m, e) / 2;
							sumb += t1.getEntry(m, a) * f.getEntry(m, e) / 2;
						}
						sum1 += t2.getEntry(i, j, a, e)
						        * (f.getEntry(b, e) - suma);
						sum2 += t2.getEntry(i, j, b, e)
						        * (f.getEntry(a, e) - sumb);
					}

					for(int m = 0; m < sofock.getShape(0); m++) {
						double suma = 0, sumb = 0;
						for(int e = occupied; e < sofock.getShape(0); e++) {
							suma += t1.getEntry(j, e) * f.getEntry(m, e) / 2;
							sumb += t1.getEntry(i, e) * f.getEntry(m, e) / 2;
						}
						sum3 += t2.getEntry(i, m, a, b)
						        * (f.getEntry(m, j) + suma);
						sum4 += t2.getEntry(j, m, a, b)
						        * (f.getEntry(m, i) + sumb);
					}

					for(int m = 0; m < sofock.getShape(0); m++) {
						for(int n = 0; n < sofock.getShape(0); n++) {
							sum5 += (t2.getEntry(m, n, a, b)
							        + t1.getEntry(m, a) * t1.getEntry(n, b)
							        - t1.getEntry(m, b) * t1.getEntry(n, a))
							        * w.getEntry(m, n, i, j) / 2;
						}
					}

					for(int e = occupied; e < sofock.getShape(0); e++) {
						for(int f = occupied; f < sofock.getShape(0); f++) {
							sum6 += (t2.getEntry(i, j, e, f)
							        + t1.getEntry(i, e) * t1.getEntry(j, f)
							        - t1.getEntry(i, f) * t1.getEntry(j, e))
							        * w.getEntry(a, b, e, f) / 2;
						}
					}

					for(int m = 0; m < occupied; m++) {
						for(int e = occupied; e < sofock.getShape(0); e++) {
							sum7 += t2.getEntry(i, m, a, e)
							        * w.getEntry(m, b, e, j)
							        - t1.getEntry(i, e) * t1.getEntry(m, a)
							                * sotei.getEntry(m, b, e, j);
							sum8 += t2.getEntry(j, m, a, e)
							        * w.getEntry(m, b, e, i)
							        - t1.getEntry(j, e) * t1.getEntry(m, a)
							                * sotei.getEntry(m, b, e, i);
							sum9 += t2.getEntry(i, m, b, e)
							        * w.getEntry(m, a, e, j)
							        - t1.getEntry(i, e) * t1.getEntry(m, b)
							                * sotei.getEntry(m, a, e, j);
							sum10 += t2.getEntry(j, m, b, e)
							        * w.getEntry(m, a, e, i)
							        - t1.getEntry(j, e) * t1.getEntry(m, b)
							                * sotei.getEntry(m, a, e, i);
						}
					}

					for(int e = occupied; e < sofock.getShape(0); e++) {
						sum11 += t1.getEntry(i, e) * sotei.getEntry(a, b, e, j);
						sum12 += t1.getEntry(j, e) * sotei.getEntry(a, b, e, i);
					}

					for(int m = 0; m < sofock.getShape(0); m++) {
						sum13 += t1.getEntry(m, a) * sotei.getEntry(m, b, i, j);
						sum14 += t1.getEntry(m, b) * sotei.getEntry(m, a, i, j);
					}

					if(sofock.getEntry(i, i) + sofock.getEntry(j, j)
					        - sofock.getEntry(a, a) - sofock.getEntry(b, b)
					        == 0) {
						throw(new std::exception());
					}

					if(fabs(
					        (sotei.getEntry(i, j, a, b) + sum1 - sum2 - sum3
					                + sum4 + sum5 + sum6 + sum7 - sum8 - sum9
					                + sum10 + sum11 - sum12 + sum13 - sum14)
					                / (sofock.getEntry(i, i)
					                        + sofock.getEntry(j, j)
					                        - sofock.getEntry(a, a)
					                        - sofock.getEntry(b, b))) >= 10) {
						throw(new std::exception());
					}

					out.setEntry(
					        (sotei.getEntry(i, j, a, b) + sum1 - sum2 - sum3
					                + sum4 + sum5 + sum6 + sum7 - sum8 - sum9
					                + sum10 + sum11 - sum12 - sum13 + sum14)
					                / (sofock.getEntry(i, i)
					                        + sofock.getEntry(j, j)
					                        - sofock.getEntry(a, a)
					                        - sofock.getEntry(b, b)), i, j, a,
					        b);
				}
			}
		}
	}
}

double compchem::strategies::DefaultCCSDCorrection::CCEnergy(
        const compchem::AbstractMatrix<double> &orbitals,
        const compchem::AbstractMatrix<double> &fock,
        const compchem::AbstractMatrix<double> &teri, int nelectrons) {
	compchem::Matrix<double> *sotei, *sofock, *f = new compchem::Matrix<double>(
	        {2 * fock.getShape(0), 2 * fock.getShape(0)}), *w =
	        new compchem::Matrix<double>(
	                {2 * fock.getShape(0), 2 * fock.getShape(0), 2
	                        * fock.getShape(0), 2 * fock.getShape(0)}), *t1 =
	        new compchem::Matrix<double>(
	                {2 * fock.getShape(0), 2 * fock.getShape(0)}), *t2 =
	        new compchem::Matrix<double>(
	                {2 * fock.getShape(0), 2 * fock.getShape(0), 2
	                        * fock.getShape(0), 2 * fock.getShape(0)}), *hold1 =
	        new compchem::Matrix<double>(
	                {2 * fock.getShape(0), 2 * fock.getShape(0)}), *hold2 =
	        new compchem::Matrix<double>(
	                {2 * fock.getShape(0), 2 * fock.getShape(0), 2
	                        * fock.getShape(0), 2 * fock.getShape(0)});

	double energy = 0, energy_prev = 1;
	sotei = &calculateSOTEI(teri, orbitals, nelectrons);
	sofock = &calculateSOFock(fock, orbitals);

	for(int i = 0; i < nelectrons; i++) {
		for(int a = nelectrons; a < sofock->getShape(0); a++) {
			t1->setEntry(0, i, a);
		}
	}

	for(int i = 0; i < nelectrons; i++) {
		for(int j = 0; j < nelectrons; j++) {
			for(int a = nelectrons; a < sofock->getShape(0); a++) {
				for(int b = nelectrons; b < sofock->getShape(0); b++) {
					t2->setEntry(
					        sotei->getEntry(i, j, a, b)
					                / (sofock->getEntry(i, i)
					                        + sofock->getEntry(j, j)
					                        - sofock->getEntry(a, a)
					                        - sofock->getEntry(b, b)), i, j, a,
					        b);
				}
			}
		}
	}

	double sum = 0;

	for(int i = 0; i < nelectrons; i++) {
		for(int j = 0; j < nelectrons; j++) {
			for(int a = nelectrons; a < sofock->getShape(0); a++) {
				for(int b = nelectrons; b < sofock->getShape(0); b++) {
					sum += sotei->getEntry(i, j, a, b)
					        * t2->getEntry(i, j, a, b);
				}
			}
		}
	}
	energy = sum / 4;

	while(fabs(energy - energy_prev) > 0.000001) {
		energy_prev = energy;
		calculateFIntermediate(*sofock, *sotei, *t1, *t2, nelectrons, *f);
		calculateWIntermediate(*sofock, *sotei, *t1, *t2, nelectrons, *w);
		calculateT1(*sofock, *sotei, *t1, *t2, *f, *w, nelectrons, *hold1);
		calculateT2(*sofock, *sotei, *t1, *t2, *f, *w, nelectrons, *hold2);

		for(int i = 0; i < sofock->getShape(0); i++) {
			for(int j = 0; j < sofock->getShape(0); j++) {
				t1->setEntry(hold1->getEntry(i, j), i, j);
				for(int k = 0; k < sofock->getShape(0); k++) {
					for(int l = 0; l < sofock->getShape(0); l++) {
						t2->setEntry(hold2->getEntry(i, j, k, l), i, j, k, l);
					}
				}
			}
		}

		energy = 0;
		double sum1 = 0, sum2 = 0, sum3 = 0;
		for(int i = 0; i < nelectrons; i++) {
			for(int a = nelectrons; a < sofock->getShape(0); a++) {
				if(std::isnan(t1->getEntry(i, a))) {
					throw(new std::exception());
				}
				sum1 += sofock->getEntry(i, a) * t1->getEntry(i, a);
			}
		}

		for(int i = 0; i < nelectrons; i++) {
			for(int j = 0; j < nelectrons; j++) {
				for(int a = nelectrons; a < sofock->getShape(0); a++) {
					for(int b = nelectrons; b < sofock->getShape(0); b++) {
						if(std::isnan(t2->getEntry(i, j, a, b))) {
							throw(new std::exception());
						}
						sum2 += sotei->getEntry(i, j, a, b)
						        * t2->getEntry(i, j, a, b) / 4;
						sum3 += sotei->getEntry(i, j, a, b) * t1->getEntry(i, a)
						        * t1->getEntry(j, b) / 2;
					}
				}
			}
		}
		energy = sum1 + sum2 + sum3;
	}

	delete f;
	delete w;
	delete hold1;
	delete hold2;
	delete t1;
	delete t2;
	delete sofock;
	delete sotei;

	return (energy);
}

