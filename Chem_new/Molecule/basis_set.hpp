/*
 * basis_set.hpp
 *
 *  Created on: Feb 11, 2019
 *      Author: connor
 */

#ifndef BASIS_SET_HPP_
#define BASIS_SET_HPP_

#include "molecule.hpp"

namespace compchem {

class AbstractBasisSet {
public:
	AbstractBasisSet() {
		;
	}
	virtual ~AbstractBasisSet() {
		;
	}

	virtual int norbitals(const compchem::AbstractMolecule &mol) const = 0;
	virtual int norbitals(int z) const = 0;
};

}



#endif /* BASIS_SET_HPP_ */
