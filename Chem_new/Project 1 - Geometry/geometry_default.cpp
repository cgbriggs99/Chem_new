/*
 * geometry_default.cpp
 *
 *  Created on: Jan 28, 2019
 *      Author: cgbri
 */

#include "geometry_default.hpp"
#include "../Base/math.hpp"

compchem::Matrix<double> &
compchem::strategies::DefaultGeometryStrategy::findDistances(compchem::AbstractMolecule &mol) {
	compchem::Matrix<double> *out =
			new compchem::Matrix<double>({mol.getNumAtoms(), mol.getNumAtoms()});

	std::vector<Atom> &vec = mol.getAtoms();
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
compchem::strategies::DefaultGeometryStrategy::findBondAngles(AbstractMolecule &mol) {
	compchem::Matrix<double> *out =
			new compchem::Matrix<double>({mol.getNumAtoms(),
		mol.getNumAtoms(), mol.getNumAtoms()});
	std::vector<Atom> &vec = mol.getAtoms();
	for(int i = 0; i < mol.getNumAtoms(); i++) {
		for(int j = 0; j < mol.getNumAtoms(); j++) {
			for(int k = 0; k < j; k++) {
				if(i == j || i == k || j == k) {
					out->getEntry({i, j, k}) = 0;
					continue;
				}
				std::vector<double> vij = {vec[i].getX() - vec[j].getX(),
				                           vec[i].getY() - vec[j].getY(),
				                           vec[i].getZ() - vec[j].getZ()
				}, vik = {vec[i].getX() - vec[k].getX(),
				          vec[i].getY() - vec[k].getY(),
				          vec[i].getZ() - vec[k].getZ()
				};
				double rij = hypot3(vij[0], vij[1], vij[2]), rik = hypot3(vik[0], vik[1], vik[2]);
				double dot = dotprod(vij, vik);
				out->getEntry({i, j, k}) = acos(dot / (rij * rik));
				out->getEntry({i, k, j}) = out->getEntry({i, k , j});
			}
		}
	}
	return (*out);
}

compchem::Matrix<double> &
compchem::strategies::DefaultGeometryStrategy::findPlaneAngles(AbstractMolecule &mol) {
	Matrix<double> *out = new Matrix<double>({mol.getNumAtoms(), mol.getNumAtoms(),
		mol.getNumAtoms(), mol.getNumAtoms()});
	std::vector<Atom> &vec = mol.getAtoms();
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
					std::vector<double> vik = {vec[i].getX() - vec[k].getX(),
					                           vec[i].getY() - vec[k].getY(),
					                           vec[i].getZ() - vec[k].getZ()
					};
					std::vector<double> vil = {vec[i].getX() - vec[l].getX(),
					                           vec[i].getY() - vec[l].getY(),
					                           vec[i].getZ() - vec[l].getZ()
					};
					double rij = hypot3(vij[0], vij[1], vij[2]);
					double rik = hypot3(vik[0], vik[1], vik[2]);
					double ril = hypot3(vil[0], vil[1], vil[2]);
					double angle = mol.getBondAngles().getEntry({i, j, k});
					out->getEntry({i, j, k, l}) = asin((dotprod(crossprod(vij, vik), vil)) / (rij * rik * ril * sin(angle)));
				}
			}
		}
	}
	return (*out);
}

compchem::Matrix<double> &
compchem::strategies::DefaultGeometryStrategy::findPlaneAngles(AbstractMolecule &mol) {
	Matrix<double> *out = new Matrix<double>({mol.getNumAtoms(), mol.getNumAtoms(),
		mol.getNumAtoms(), mol.getNumAtoms()});
	std::vector<Atom> &vec = mol.getAtoms();
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
					double angle1 = mol.getBondAngles().getEntry({j, i, k});
					double angle2 = mol.getBondAngles().getEntry({k, j, l});
					out->getEntry({i, j, k, l}) = acos((dotprod(crossprod(vij, vjk), crossprod(vjk, vkl))) /
							(rij * rjk * rjk * rkl * sin(angle1) * sin(angle2)));
				}
			}
		}
	}
	return (*out);
}

std::vector<double> &
compchem::strategies::DefaultGeometryStrategy::findCenterOfMass(compchem::AbstractMolecule &mol) {
	std::vector<double> *out = new std::vector<double>({0, 0, 0});

	std::vector<Atom> &vec = mol.getAtoms();
	double mass = 0;
	for(int i = 0; i < mol.getNumAtoms(); i++) {
		mass += vec[i].getMass();
		out->at(0) += vec[i].pos[0];
		out->at(1) += vec[i].pos[1];
		out->at(2) += vec[i].pos[2];
	}
	out->at(0) /= mass;
	out->at(1) /= mass;
	out->at(2) /= mass;

	return (*out);
}


std::vector<double> &
compchem::strategies::DefaultGeometryStrategy::findRotationalConstants(compchem::AbstractMolecule &mol) {
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
	std::vector<Atom> &vec = mol.getAtoms();
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

	Matrix<double> &temp = eigenstrat->eigenvals(*moments);
	std::vector<double> out;
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

	delete moments;
	delete &temp;

	return (out);
}

compchem::rotor_type compchem::strategies::DefaultGeometryStrategy::findRotor(compchem::AbstractMolecule &mol) {
	std::vector<double> &moms = mol.getMoments();
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


