/*
 * geometry_default.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: cgbri
 */

#ifndef GEOMETRY_DEFAULT_HPP_
#define GEOMETRY_DEFAULT_HPP_

#include "geometry.hpp"
#include "../Base/eigenvalues.hpp"

namespace compchem {

namespace strategies {

class DefaultGeometryStrategy : public GeometryCalcStrategy {
private:
	EigenvalueStrategy<double> *eigenstrat;
public:
	DefaultGeometryStrategy() {
		;
	}

	~DefaultGeometryStrategy() {
		;
	}

	compchem::Matrix<double> &findDistances(compchem::AbstractMolecule &mol) override;
	compchem::Matrix<double> &findBondAngles(compchem::AbstractMolecule &mol) override;
	compchem::Matrix<double> &findPlaneAngles(compchem::AbstractMolecule &mol) override;
	compchem::Matrix<double> &findTorsionAngles(compchem::AbstractMolecule &mol) override;
	std::vector<double> &findRotationalConstants(compchem::AbstractMolecule &mol) override;
	std::vector<double> &findCenterOfMass(compchem::AbstractMolecule &mol) override;
	compchem::rotor_type findRotor(compchem::AbstractMolecule &mol) override;

};

}
}

#endif /* GEOMETRY_DEFAULT_HPP_ */
