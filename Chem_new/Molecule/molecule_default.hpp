/*
 * molecule_default.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: connor
 */

#ifndef MOLECULE_DEFAULT_HPP_
#define MOLECULE_DEFAULT_HPP_

#include "molecule.hpp"
#include "../Base/matrix_default.hpp"

namespace compchem {
namespace strategies {

class DefaultMolecule : public AbstractMolecule {
private:
	std::vector<Atom> *atoms;
	int num_electrons, num_orbitals;
public:
	DefaultMolecule() {
		atoms = new std::vector<Atom>();
		num_electrons = 0;
		num_orbitals = 0;
	}

	//Defined in molecule.cpp
	~DefaultMolecule();

	const std::vector<Atom> &getAtoms() const override {
		return (*(this->atoms));
	}

	void addAtom(Atom a) {
		this->atoms->push_back(a);
		num_electrons += a.getAtomicNum();
		num_orbitals += compchem::orbitals((int) a.getAtomicNum());
	}

	int norbital() const {
		return (num_orbitals);
	}

	int nelectron() const {
		return (num_electrons);
	}

	void translateAtoms(const std::vector<double> &diff) override {
		for(int i = 0; i < this->natom(); i++) {
			this->atoms->at(i).setX(this->atoms->at(i).getX() + diff[0]);
			this->atoms->at(i).setY(this->atoms->at(i).getY() + diff[1]);
			this->atoms->at(i).setZ(this->atoms->at(i).getZ() + diff[2]);
		}
	}

	void translateCOM(const std::vector<double> &diff) override {
		for(int i = 0; i < this->natom(); i++) {
			this->atoms->at(i).setX(this->atoms->at(i).getX() - diff[0]);
			this->atoms->at(i).setY(this->atoms->at(i).getY() - diff[1]);
			this->atoms->at(i).setZ(this->atoms->at(i).getZ() - diff[2]);
		}
	}

	double x(int i) const override {
		return (this->atoms->at(i).getX());
	}

	double y(int i) const override {
		return (this->atoms->at(i).getY());
	}

	double z(int i) const override {
		return (this->atoms->at(i).getZ());
	}

	double mass(int i) const override {
		return (this->atoms->at(i).getMass());
	}

	double fZ(int i) const override {
		return (this->atoms->at(i).getAtomicNum());
	}
};

}
}



#endif /* MOLECULE_DEFAULT_HPP_ */
