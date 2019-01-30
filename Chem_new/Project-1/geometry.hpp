/*
 * geometry.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: cgbri
 */

#ifndef GEOMETRY_HPP_
#define GEOMETRY_HPP_

#include "../Base/matrix.hpp"
//#include "../Molecule/molecule.hpp"
//#include <psi4/libmints/molecule.h>
#include "/usr/local/psi4/include/psi4/pragma.h"
#include "/usr/local/psi4/include/psi4/libmints/molecule.h"

namespace compchem {


class GeometryCalcStrategy {
public:
	GeometryCalcStrategy() = default;
	virtual ~GeometryCalcStrategy() = default;

	virtual compchem::Matrix<double> &findDistances(const psi::Molecule &mol) = 0;
	virtual compchem::Matrix<double> &findBondAngles(const psi::Molecule &mol) = 0;
	virtual compchem::Matrix<double> &findPlaneAngles(const psi::Molecule &mol, const compchem::Matrix<double> &bond_angles) = 0;
	virtual compchem::Matrix<double> &findTorsionAngles(const psi::Molecule &mol, const compchem::Matrix<double> &bond_angles) = 0;
	virtual std::vector<double> &findCenterOfMass(const psi::Molecule &mol) = 0;
	virtual std::vector<double> &findPrincipleMoments(const compchem::Matrix<double> &moment) = 0;
//	virtual std::vector<double> &findRotationalConstants(const compchem::AbstractMolecule &mol) = 0;
	virtual compchem::Matrix<double> &findMoments(const psi::Molecule &mol) = 0;
	virtual psi::RotorType findRotor(const std::vector<double> &principles) = 0;
};

}

#endif /* GEOMETRY_HPP_ */
