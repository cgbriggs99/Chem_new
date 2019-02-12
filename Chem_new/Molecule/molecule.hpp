/*
 * molecule.hpp
 *
 *  Created on: Dec 22, 2018
 *      Author: connor
 */

#ifndef MOLECULE_HPP_
#define MOLECULE_HPP_

#include <vector>
#include "../Base/matrix.hpp"

namespace compchem {

double amu(int z);
int orbitals(int n);

typedef enum {
	SPHERICAL, LINEAR, OBLATE, PROLATE, ASYMMETRIC
} rotor_type;

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

	double getX() const {
		return (pos[0]);
	}

	double getY() const {
		return (pos[1]);
	}

	double getZ() const {
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
		return (atomic_num);
	}

	void setAtomicNum(int atomicNum) {
		atomic_num = atomicNum;
	}

	int getCharge() const {
		return (charge);
	}

	void setCharge(int charge) {
		this->charge = charge;
	}

	double getMass() const {
		return (mass);
	}

	void setMass(double mass) {
		this->mass = mass;
	}

	double getTrueCharge() const {
		return (true_charge);
	}

	void setTrueCharge(double trueCharge) {
		true_charge = trueCharge;
	}
};

class AbstractMolecule {
public:
	AbstractMolecule() {
		;
	}
	virtual ~AbstractMolecule() = default;

	/*
	 * My body for the molecule class. Will become deprecated.
	 */
	virtual const std::vector<Atom> &getAtoms() const = 0;
	virtual int natom() const {
		return (this->getAtoms().size());
	}
	virtual int nelectron() const = 0;
	virtual void addAtom(Atom a) = 0;

	virtual void translateAtoms(const std::vector<double> &diff) = 0;
	virtual void translateCOM(const std::vector<double> &diff) = 0;

	virtual double x(int i) const = 0;
	virtual double y(int i) const = 0;
	virtual double z(int i) const = 0;
	virtual double mass(int i) const = 0;
	virtual double fZ(int i) const = 0;
};
}

#endif /* MOLECULE_HPP_ */
