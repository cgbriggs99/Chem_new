/*
 * sto3g_basis_set.hpp
 *
 *  Created on: Feb 11, 2019
 *      Author: connor
 */

#ifndef STO3G_BASIS_SET_HPP_
#define STO3G_BASIS_SET_HPP_

#include "basis_set.hpp"

namespace compchem {
namespace strategies {

class STO3GBasisSet : public AbstractBasisSet {
public:
	STO3GBasisSet() {
		;
	}
	virtual ~STO3GBasisSet() {
		;
	}

	int norbitals(const compchem::AbstractMolecule &mol) const override {
		int sum = 0;
		for(int i = 0; i < mol.natom(); i++) {
			sum += compchem::orbitals(mol.fZ(i));
		}
		return (sum);
	}

	int norbitals(int z) const override {
		return (compchem::orbitals(z));
	}
};

}
}



#endif /* STO3G_BASIS_SET_HPP_ */
