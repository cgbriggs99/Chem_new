/*
 * molecule.cpp
 *
 *  Created on: Dec 22, 2018
 *      Author: connor
 */

#include "molecule_default.hpp"

compchem::strategies::DefaultMolecule::~DefaultMolecule()  {
	if(dists != nullptr) {
		delete dists;
	}
	if(bonds != nullptr) {
		delete bonds;
	}
	if(plane_angles != nullptr) {
		delete plane_angles;
	}
	if(torsion != nullptr) {
		delete torsion;
	}
}

