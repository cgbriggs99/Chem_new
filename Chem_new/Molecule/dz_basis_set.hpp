/*
 * dz_basis_set.hpp
 *
 *  Created on: Feb 11, 2019
 *      Author: connor
 */

#ifndef DZ_BASIS_SET_HPP_
#define DZ_BASIS_SET_HPP_

#include "basis_set.hpp"

namespace compchem {
namespace strategies {

class DZBasisSet : public AbstractBasisSet {
public:
	DZBasisSet() {
		;
	}
	virtual ~DZBasisSet() {
		;
	}

	int norbitals(const compchem::AbstractMolecule &mol) const override {
		int sum = 0;
		for(int i = 0; i < mol.natom(); i++) {
			sum += compchem::orbitals(mol.fZ(i));
		}
		return (2 * sum);
	}

	int norbitals(int z) const override {
		return (2 * compchem::orbitals(z));
	}
};

}
}



#endif /* DZ_BASIS_SET_HPP_ */
