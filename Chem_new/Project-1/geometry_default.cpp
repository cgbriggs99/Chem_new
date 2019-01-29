/*
 * geometry_default.cpp
 *
 *  Created on: Jan 28, 2019
 *      Author: cgbri
 */

#include "../Project-1/geometry_default.hpp"

#include "../Base/math.hpp"

compchem::Matrix<double> &
compchem::strategies::DefaultGeometryStrategy::findDistances(const compchem::AbstractMolecule &mol) {
	compchem::Matrix<double> *out =
			new compchem::Matrix<double>({mol.getNumAtoms(), mol.getNumAtoms()});

	const std::vector<Atom> &vec = mol.getAtoms();
	for(int i = 0; i < mol.getNumAtoms(); i++) {
		out->getEntry({i, i}) = 0;
		for(int j = 0; j < i; j++) {
			double dx = vec[i].getX() - vec[j].getX();
			double dy = vec[i].getY() - vec[j].getY();
			double dz = vec[i].getZ() - vec[j].getZ();
			out->getEntry({i, j}) = hypot3(dx, dy, dz);
			out->getEntry({j, i}) = out->getEntry({i, j});
		}
	}
	return (*out);
}

compchem::Matrix<double> &
compchem::strategies::DefaultGeometryStrategy::findBondAngles(const AbstractMolecule &mol) {
	compchem::Matrix<double> *out =
			new compchem::Matrix<double>({mol.getNumAtoms(),
		mol.getNumAtoms(), mol.getNumAtoms()});
	const std::vector<Atom> &vec = mol.getAtoms();
	for(int i = 0; i < mol.getNumAtoms(); i++) {
		for(int j = 0; j < mol.getNumAtoms(); j++) {
			for(int k = 0; k < mol.getNumAtoms(); k++) {
				if(i == j || i == k || j == k) {
					out->getEntry({i, j, k}) = 0;
					continue;
				}
				std::vector<double> vji = {vec[j].getX() - vec[i].getX(),
				                           vec[j].getY() - vec[i].getY(),
				                           vec[j].getZ() - vec[i].getZ()
				}, vjk = {vec[j].getX() - vec[k].getX(),
				          vec[j].getY() - vec[k].getY(),
				          vec[j].getZ() - vec[k].getZ()
				};
				double rij = hypot3(vji[0], vji[1], vji[2]), rjk = hypot3(vjk[0], vjk[1], vjk[2]);
				double dot = dotprod(vji, vjk);
				out->getEntry({i, j, k}) = acos(dot / (rij * rjk));
			}
		}
	}
	return (*out);
}

compchem::Matrix<double> &
compchem::strategies::DefaultGeometryStrategy::findPlaneAngles(const AbstractMolecule &mol) {
	Matrix<double> *out = new Matrix<double>({mol.getNumAtoms(), mol.getNumAtoms(),
		mol.getNumAtoms(), mol.getNumAtoms()});
	const std::vector<Atom> &vec = mol.getAtoms();
	for(int i = 0; i < mol.getNumAtoms(); i++) {
		for(int j = 0; j < mol.getNumAtoms(); j++) {
			for(int k = 0; k < mol.getNumAtoms(); k++) {
				for(int l = 0; l < mol.getNumAtoms(); l++) {
					if(i == j || i == k || i == l || j == k || j == l || k == l) {
						out->getEntry({i, j, k, l}) = 0;
						continue;
					}
					std::vector<double> vkj = {vec[k].getX() - vec[j].getX(),
					                           vec[k].getY() - vec[j].getY(),
					                           vec[k].getZ() - vec[j].getZ()
					};
					std::vector<double> vki = {vec[k].getX() - vec[i].getX(),
					                           vec[k].getY() - vec[i].getY(),
					                           vec[k].getZ() - vec[i].getZ()
					};
					std::vector<double> vkl = {vec[k].getX() - vec[l].getX(),
					                           vec[k].getY() - vec[l].getY(),
					                           vec[k].getZ() - vec[l].getZ()
					};
					double rkj = hypot3(vkj[0], vkj[1], vkj[2]);
					double rki = hypot3(vki[0], vki[1], vki[2]);
					double rkl = hypot3(vkl[0], vkl[1], vkl[2]);

					std::vector<double> ekj = {vkj[0] / rkj, vkj[1] / rkj, vkj[2] / rkj};
					std::vector<double> eki = {vki[0] / rki, vki[1] / rki, vki[2] / rki};
					std::vector<double> ekl = {vkl[0] / rkl, vkl[1] / rkl, vkl[2] / rkl};
					double angle = mol.getBondAngles().getEntry({j, k, l});
					double hold = (dotprod(crossprod(ekl, ekj), eki)) / (sin(angle));
					if(hold >= 1) {
						out->getEntry({i, j, k, l}) = M_PI_2;
					} else if(hold <= -1) {
						out->getEntry({i, j, k, l}) = -M_PI_2;
					} else {
						out->getEntry({i, j, k, l}) = asin(hold);
					}
				}
			}
		}
	}
	return (*out);
}

compchem::Matrix<double> &
compchem::strategies::DefaultGeometryStrategy::findTorsionAngles(const compchem::AbstractMolecule &mol) {
	Matrix<double> *out = new Matrix<double>({mol.getNumAtoms(), mol.getNumAtoms(),
		mol.getNumAtoms(), mol.getNumAtoms()});
	const std::vector<Atom> &vec = mol.getAtoms();
	for(int i = 0; i < mol.getNumAtoms(); i++) {
		for(int j = 0; j < mol.getNumAtoms(); j++) {
			for(int k = 0; k < mol.getNumAtoms(); k++) {
				for(int l = 0; l < mol.getNumAtoms(); l++) {
					if(i == j || i == k || i == l || j == k || j == l || k == l) {
						out->getEntry({i, j, k, l}) = 0;
						continue;
					}
					std::vector<double> vij = {vec[i].getX() - vec[j].getX(),
					                           vec[i].getY() - vec[j].getY(),
					                           vec[i].getZ() - vec[j].getZ()
					};
					std::vector<double> vjk = {vec[j].getX() - vec[k].getX(),
					                           vec[j].getY() - vec[k].getY(),
					                           vec[j].getZ() - vec[k].getZ()
					};
					std::vector<double> vkl = {vec[k].getX() - vec[l].getX(),
					                           vec[k].getY() - vec[l].getY(),
					                           vec[k].getZ() - vec[l].getZ()
					};
					double rij = hypot3(vij[0], vij[1], vij[2]);
					double rjk = hypot3(vjk[0], vjk[1], vjk[2]);
					double rkl = hypot3(vkl[0], vkl[1], vkl[2]);

//					std::vector<double> eij = {vij[0] / rij, vij[1] / rij, vij[2] / rij};
//					std::vector<double> ejk = {vjk[0] / rjk, vjk[1] / rjk, vjk[2] / rjk};
//					std::vector<double> ekl = {vkl[0] / rkl, vkl[1] / rkl, vkl[2] / rkl};

					double angle1 = mol.getBondAngles().getEntry({i, j, k});
					double angle2 = mol.getBondAngles().getEntry({j, k, l});
					double hold = (dotprod(crossprod(vjk, vij), crossprod(vkl, vjk))) /
							(rij * rjk * rjk * rkl * sin(angle1) * sin(angle2));
					if(hold >= 1) {
						out->getEntry({i, j, k, l}) = 0;
					} else if(hold <= -1) {
						out->getEntry({i, j, k, l}) = M_PI;
					} else {
						out->getEntry({i, j, k, l}) = acos(hold);
					}
				}
			}
		}
	}
	return (*out);
}

std::vector<double> &
compchem::strategies::DefaultGeometryStrategy::findCenterOfMass(const compchem::AbstractMolecule &mol) {
	std::vector<double> *out = new std::vector<double>({0, 0, 0});

	const std::vector<Atom> &vec = mol.getAtoms();
	double mass = 0;
	for(int i = 0; i < mol.getNumAtoms(); i++) {
		mass += vec[i].getMass();
		out->at(0) += vec[i].getMass() * vec[i].getX();
		out->at(1) += vec[i].getMass() * vec[i].getY();
		out->at(2) += vec[i].getMass() * vec[i].getZ();
	}
	out->at(0) /= mass;
	out->at(1) /= mass;
	out->at(2) /= mass;

	return (*out);
}


compchem::Matrix<double> &compchem::strategies::DefaultGeometryStrategy::findMoments(
		const compchem::AbstractMolecule &mol) {
	Matrix<double> *moments = new Matrix<double>({3, 3});
	moments->getEntry({0, 0}) = 0;
	moments->getEntry({1, 0}) = 0;
	moments->getEntry({2, 0}) = 0;
	moments->getEntry({0, 1}) = 0;
	moments->getEntry({1, 1}) = 0;
	moments->getEntry({2, 1}) = 0;
	moments->getEntry({0, 2}) = 0;
	moments->getEntry({1, 2}) = 0;
	moments->getEntry({2, 2}) = 0;
	const std::vector<Atom> &vec = mol.getAtoms();
	for(int i = 0; i < mol.getNumAtoms(); i++) {
		double x = vec[i].getX(), y = vec[i].getY(), z = vec[i].getZ(), mass = vec[i].getMass();
		moments->getEntry({0, 0}) += mass * (y * y + z * z);
		moments->getEntry({1, 0}) += mass * x * y;
		moments->getEntry({2, 0}) += mass * x * z;
		moments->getEntry({0, 1}) += mass * x * y;
		moments->getEntry({1, 1}) += mass * (x * x + z * z);
		moments->getEntry({2, 1}) += mass * y * z;
		moments->getEntry({0, 2}) += mass * x * z;
		moments->getEntry({1, 2}) += mass * y * z;
		moments->getEntry({2, 2}) += mass * (x * x + y * y);
	}
	return (*moments);
}
std::vector<double> &
compchem::strategies::DefaultGeometryStrategy::findPrincipleMoments(const compchem::AbstractMolecule &mol) {
	Matrix<double> &temp = eigenstrat->eigenvals(mol.getMoments());
	std::vector<double> &out = *new std::vector<double>(3);
	double a, b, c;
	a = temp.getEntry({0});
	b = temp.getEntry({1});
	c = temp.getEntry({2});

	if(a <= b && a <= c) {
		out[0] = a;
		if(b <= c) {
			out[1] = b;
			out[2] = c;
		} else {
			out[1] = c;
			out[2] = b;
		}
	} else if(b <= a && b <= c) {
		out[0] = b;
		if(a <= c) {
			out[1] = a;
			out[2] = c;
		} else {
			out[1] = c;
			out[2] = a;
		}
	} else {
		out[0] = c;
		if(a <= b) {
			out[1] = a;
			out[2] = b;
		} else {
			out[1] = b;
			out[2] = a;
		}
	}

	delete &temp;

	return (out);
}

compchem::rotor_type compchem::strategies::DefaultGeometryStrategy::findRotor(const compchem::AbstractMolecule &mol) {
	const std::vector<double> &moms = mol.getPrincipleMoments();
	if(compareDoubles(moms[0], moms[1], 0.0001) == 0 && compareDoubles(moms[1], moms[2], 0.0001) == 0) {
		return (SPHERICAL);
	}
	if(compareDoubles(moms[0], 0, 0.0001) == 0) {
		return (LINEAR);
	}
	if(compareDoubles(moms[0], moms[1], 0.0001) == 0) {
		return (OBLATE);
	}
	if(compareDoubles(moms[1], moms[2], 0.0001) == 0) {
		return (PROLATE);
	}
	return (ASYMMETRIC);
}


