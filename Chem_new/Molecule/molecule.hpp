/*
 * molecule.hpp
 *
 *  Created on: Dec 22, 2018
 *      Author: connor
 */

#ifndef MOLECULE_HPP_
#define MOLECULE_HPP_

#include <vector>

namespace compchem {

class Atom {
private:
	double pos[3];
	double mass;
	double true_charge;
	int charge;
	int atomic_num;
public:
	Atom(int number, double mass, int charge, double x, double y, double z) {
		atomic_num = number;
		this->mass = mass;
		pos[0] = x;
		pos[1] = y;
		pos[2] = z;
		this->charge = charge;
		true_charge = charge;
	}

	double getX() {
		return (pos[0]);
	}

	double getY() {
		return (pos[1]);
	}

	double getZ() {
		return (pos[2]);
	}

	void setX(double x) {
		pos[0] = x;
	}

	void setY(double y) {
		pos[1] = y;
	}

	void setZ(double z) {
		pos[2] = z;
	}

	int getAtomicNum() const {
		return atomic_num;
	}

	void setAtomicNum(int atomicNum) {
		atomic_num = atomicNum;
	}

	int getCharge() const {
		return charge;
	}

	void setCharge(int charge) {
		this->charge = charge;
	}

	double getMass() const {
		return mass;
	}

	void setMass(double mass) {
		this->mass = mass;
	}

	double getTrueCharge() const {
		return true_charge;
	}

	void setTrueCharge(double trueCharge) {
		true_charge = trueCharge;
	}
};

class AbstractMolecule {
public:
	virtual ~AbstractMolecule();

	virtual std::vector<Atom> getAtoms();
	virtual void addAtom(Atom a);
};
}

#endif /* MOLECULE_HPP_ */
