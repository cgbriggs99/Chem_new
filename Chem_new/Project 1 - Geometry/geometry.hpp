/*
 * geometry.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: cgbri
 */

#ifndef GEOMETRY_HPP_
#define GEOMETRY_HPP_

#include "../Base/matrix.hpp"
#include "../Molecule/molecule.hpp"

namespace compchem {


class GeometryCalcStrategy {
public:
	GeometryCalcStrategy() = default;
	virtual ~GeometryCalcStrategy() = default;

	virtual compchem::Matrix<double> &findDistances(compchem::AbstractMolecule &mol) = 0;
	virtual compchem::Matrix<double> &findBondAngles(compchem::AbstractMolecule &mol) = 0;
	virtual compchem::Matrix<double> &findPlaneAngles(compchem::AbstractMolecule &mol) = 0;
	virtual compchem::Matrix<double> &findTorsionAngles(compchem::AbstractMolecule &mol) = 0;
	virtual std::vector<double> &findRotationalConstants(compchem::AbstractMolecule &mol) = 0;
	virtual compchem::rotor_type findRotor(compchem::AbstractMolecule &mol) = 0;
};

}

#endif /* GEOMETRY_HPP_ */
