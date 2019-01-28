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
				}
			}
		}
	}
}
