/*
 * geometry_default.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: cgbri
 */

#ifndef GEOMETRY_DEFAULT_HPP_
#define GEOMETRY_DEFAULT_HPP_

#include "/usr/local/psi4/include/psi4/pragma.h"
#include "../Base/eigenvalues.hpp"
#include "../Base/eigenvalues_default.hpp"
#include "../Project-1/geometry.hpp"

#include "/usr/local/psi4/include/psi4/libmints/molecule.h"

namespace compchem {

namespace strategies {

class DefaultGeometryStrategy : public GeometryCalcStrategy {
private:
	EigenvalueStrategy<double> *eigenstrat;
public:
	DefaultGeometryStrategy() {
		this->eigenstrat = new compchem::strategies::LapackEigenvalues<double>();
	}

	~DefaultGeometryStrategy() {
		delete eigenstrat;
	}

	compchem::Matrix<double> &findDistances(const psi::Molecule &mol) override;
	compchem::Matrix<double> &findBondAngles(const psi::Molecule &mol) override;
	compchem::Matrix<double> &findPlaneAngles(const psi::Molecule &mol, const compchem::Matrix<double> &bond_angles) override;
	compchem::Matrix<double> &findTorsionAngles(const psi::Molecule &mol, const compchem::Matrix<double> &bond_angles) override;
	std::vector<double> &findPrincipleMoments(const compchem::Matrix<double> &moment) override;
//	std::vector<double> &findRotationalConstants(const psi::Molecule &mol) override;
	std::vector<double> &findCenterOfMass(const psi::Molecule &mol) override;
	psi::RotorType findRotor(const std::vector<double> &principles) override;
	compchem::Matrix<double> &findMoments(const psi::Molecule &mol) override;

};

}
}

#endif /* GEOMETRY_DEFAULT_HPP_ */
