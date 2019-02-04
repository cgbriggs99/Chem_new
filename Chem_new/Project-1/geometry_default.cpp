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
			new compchem::Matrix<double>({mol.natom(), mol.natom()});

	for(int i = 0; i < mol.natom(); i++) {
		out->getEntry({i, i}) = 0;
		for(int j = 0; j < i; j++) {
			double dx = mol.x(i) - mol.x(j);
			double dy = mol.y(i) - mol.y(j);
			double dz = mol.z(i) - mol.z(j);
			out->getEntry({i, j}) = hypot3(dx, dy, dz);
			out->getEntry({j, i}) = out->getEntry({i, j});
		}
	}
	return (*out);
}

compchem::Matrix<double> &
compchem::strategies::DefaultGeometryStrategy::findBondAngles(const compchem::AbstractMolecule &mol) {
	compchem::Matrix<double> *out =
			new compchem::Matrix<double>({mol.natom(),
		mol.natom(), mol.natom()});
	for(int i = 0; i < mol.natom(); i++) {
		for(int j = 0; j < mol.natom(); j++) {
			for(int k = 0; k < mol.natom(); k++) {
				if(i == j || i == k || j == k) {
					out->getEntry({i, j, k}) = 0;
					continue;
				}
				double vji[] = {mol.x(j) - mol.x(i),
						mol.y(j) - mol.y(i),
						mol.z(j) - mol.z(i)
				};
				double vjk[] = {mol.x(j) - mol.x(k),
				          mol.y(j) - mol.y(k),
				          mol.z(j) - mol.z(k)
				};
				double rij = hypot3(vji[0], vji[1], vji[2]), rjk = hypot3(vjk[0], vjk[1], vjk[2]);
				double dot = vji[0] * vjk[0] + vji[1] * vjk[1] + vji[2] * vjk[2];
				out->getEntry({i, j, k}) = acos(dot / (rij * rjk));
			}
		}
	}
	return (*out);
}

compchem::Matrix<double> &
compchem::strategies::DefaultGeometryStrategy::findPlaneAngles(const compchem::AbstractMolecule &mol, const compchem::Matrix<double> &bond_angles) {
	Matrix<double> *out = new Matrix<double>({mol.natom(), mol.natom(),
		mol.natom(), mol.natom()});

	for(int i = 0; i < mol.natom(); i++) {
		for(int j = 0; j < mol.natom(); j++) {
			for(int k = 0; k < mol.natom(); k++) {
				for(int l = 0; l < mol.natom(); l++) {
					if(i == j || i == k || i == l || j == k || j == l || k == l) {
						out->getEntry({i, j, k, l}) = 0;
						continue;
					}
					double vkj[] = {mol.x(k) - mol.x(j),
					                           mol.y(k) - mol.y(j),
					                           mol.z(k) - mol.z(j)
					};
					double vki[] = {mol.x(k) - mol.x(i),
					                           mol.y(k) - mol.y(i),
					                           mol.z(k) - mol.z(i)
					};
					double vkl[] = {mol.x(k) - mol.x(l),
					                           mol.y(k) - mol.y(l),
					                           mol.z(k) - mol.z(l)
					};
					double rkj = hypot3(vkj[0], vkj[1], vkj[2]);
					double rki = hypot3(vki[0], vki[1], vki[2]);
					double rkl = hypot3(vkl[0], vkl[1], vkl[2]);

					double angle = bond_angles.getEntry({j, k, l});
					double hold = (vki[0] * (vkl[1] * vkj[2] - vkl[2] * vkj[1]) +
							vki[1] * (vkl[2] * vkj[0] - vkl[0] * vkj[2]) +
							vki[2] * (vkl[0] * vkj[1] - vkl[1] * vkj[0])) / (rkj * rki * rkl * sin(angle));
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
compchem::strategies::DefaultGeometryStrategy::findTorsionAngles(const compchem::AbstractMolecule &mol,
		const compchem::Matrix<double> &bond_angles) {
	Matrix<double> *out = new Matrix<double>({mol.natom(), mol.natom(),
		mol.natom(), mol.natom()});

	for(int i = 0; i < mol.natom(); i++) {
		for(int j = 0; j < mol.natom(); j++) {
			for(int k = 0; k < mol.natom(); k++) {
				for(int l = 0; l < mol.natom(); l++) {
					if(i == j || i == k || i == l || j == k || j == l || k == l) {
						out->getEntry({i, j, k, l}) = 0;
						continue;
					}
					double vij[] = {mol.x(i) - mol.x(j),
					                           mol.y(i) - mol.y(j),
					                           mol.z(i) - mol.z(j)
					};
					double vjk[] = {mol.x(j) - mol.x(k),
					                           mol.y(j) - mol.y(k),
					                           mol.z(j) - mol.z(k)
					};
					double vkl[] = {mol.x(k) - mol.x(l),
					                           mol.y(k) - mol.y(l),
					                           mol.z(k) - mol.z(l)
					};
					double rij = hypot3(vij[0], vij[1], vij[2]);
					double rjk = hypot3(vjk[0], vjk[1], vjk[2]);
					double rkl = hypot3(vkl[0], vkl[1], vkl[2]);

					double angle1 = bond_angles.getEntry({i, j, k});
					double angle2 = bond_angles.getEntry({j, k, l});
					double hold = ((vij[1] * vjk[2] - vij[2] * vjk[1]) * (vjk[1] * vkl[2] - vjk[2] * vkl[1]) +
							(vij[2] * vjk[0] - vij[0] * vjk[2]) * (vjk[2] * vkl[0] - vjk[0] * vkl[2]) +
							(vij[0] * vjk[1] - vij[1] * vjk[0]) * (vjk[0] * vkl[1] - vjk[1] * vkl[0])) /
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

	double mass = 0;
	for(int i = 0; i < mol.natom(); i++) {
		mass += mol.mass(i);
		out->at(0) += mol.mass(i) * mol.x(i);
		out->at(1) += mol.mass(i) * mol.y(i);
		out->at(2) += mol.mass(i) * mol.z(i);
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
	for(int i = 0; i < mol.natom(); i++) {
		double x = mol.x(i), y = mol.y(i), z = mol.z(i), mass = mol.mass(i);
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
compchem::strategies::DefaultGeometryStrategy::findPrincipleMoments(const compchem::Matrix<double> &moments) {
	Matrix<double> &temp = eigenstrat->eigenvals(moments);
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

compchem::rotor_type compchem::strategies::DefaultGeometryStrategy::findRotor(const std::vector<double> &moms) {
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


