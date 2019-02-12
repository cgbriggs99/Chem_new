/*
 * mp2.hpp
 *
 *  Created on: Feb 11, 2019
 *      Author: connor
 */

#ifndef MP2_HPP_
#define MP2_HPP_

#include "../Base/base.hpp"
#include "../Molecule/molecule.hpp"
#include "../Molecule/basis_set.hpp"

namespace compchem {

template<typename _T>
class AbstractMP2Strategy {
public:
	AbstractMP2Strategy() {
		;
	}
	virtual ~AbstractMP2Strategy() {
		;
	}

	virtual double mp2Energy(const compchem::AbstractMolecule &mol, const compchem::AbstractMatrix<double> &tei,
			const compchem::AbstractMatrix<double> &orbs, const compchem::AbstractMatrix<double> &eigs) = 0;
};

}



#endif /* MP2_HPP_ */
