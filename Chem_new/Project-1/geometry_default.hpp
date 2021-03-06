/*
 * geometry_default.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: cgbri
 */

#ifndef GEOMETRY_DEFAULT_HPP_
#define GEOMETRY_DEFAULT_HPP_

#include "../Base/eigenvalues.hpp"
#include "../Base/eigenvalues_default.hpp"
#include "../Project-1/geometry.hpp"

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

	compchem::Matrix<double> &findDistances(const compchem::AbstractMolecule &mol) override;
	compchem::Matrix<double> &findBondAngles(const compchem::AbstractMolecule &mol) override;
	compchem::Matrix<double> &findPlaneAngles(const compchem::AbstractMolecule &mol, const compchem::Matrix<double> &bond_angles) override;
	compchem::Matrix<double> &findTorsionAngles(const compchem::AbstractMolecule &mol, const compchem::Matrix<double> &bond_angles) override;
	std::vector<double> &findPrincipleMoments(const compchem::Matrix<double> &moment) override;
//	std::vector<double> &findRotationalConstants(const compchem::AbstractMolecule &mol) override;
	std::vector<double> &findCenterOfMass(const compchem::AbstractMolecule &mol) override;
	compchem::rotor_type findRotor(const std::vector<double> &principles) override;
	compchem::Matrix<double> &findMoments(const compchem::AbstractMolecule &mol) override;

};

}
}

#endif /* GEOMETRY_DEFAULT_HPP_ */
