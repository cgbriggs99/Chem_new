/*
 * geometry.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: cgbri
 */

#ifndef GEOMETRY_HPP_
#define GEOMETRY_HPP_

#include "../Base/matrix.hpp"
#include "../Base/matrix_default.hpp"
#include "../Molecule/molecule.hpp"

namespace compchem {

class GeometryCalcStrategy {
public:
	GeometryCalcStrategy() = default;
	virtual ~GeometryCalcStrategy() = default;

	virtual compchem::Matrix<double> &findDistances(
			const compchem::AbstractMolecule &mol) = 0;
	virtual compchem::Matrix<double> &findBondAngles(
			const compchem::AbstractMolecule &mol) = 0;
	virtual compchem::Matrix<double> &findPlaneAngles(
			const compchem::AbstractMolecule &mol,
			const compchem::Matrix<double> &bond_angles) = 0;
	virtual compchem::Matrix<double> &findTorsionAngles(
			const compchem::AbstractMolecule &mol,
			const compchem::Matrix<double> &bond_angles) = 0;
	virtual std::vector<double> &findCenterOfMass(
			const compchem::AbstractMolecule &mol) = 0;
	virtual std::vector<double> &findPrincipleMoments(
			const compchem::Matrix<double> &moment) = 0;
//	virtual std::vector<double> &findRotationalConstants(const compchem::AbstractMolecule &mol) = 0;
	virtual compchem::Matrix<double> &findMoments(
			const compchem::AbstractMolecule &mol) = 0;
	virtual compchem::rotor_type findRotor(
			const std::vector<double> &principles) = 0;
};

}

#endif /* GEOMETRY_HPP_ */
