/*
 * molecule_default.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: connor
 */

#ifndef MOLECULE_DEFAULT_HPP_
#define MOLECULE_DEFAULT_HPP_

#include "molecule.hpp"

namespace compchem {
namespace strategies {

class DefaultMolecule : public AbstractMolecule {
private:
	std::vector<Atom> *atoms;

	Matrix<double> *dists;
	Matrix<double> *bonds;
	Matrix<double> *plane_angles;
	Matrix<double> *torsion;
	Matrix<double> *moments;
	std::vector<double> *principle;
	std::vector<double> *rotations;
	rotor_type rotor;
public:
	DefaultMolecule() {
		dists = nullptr;
		bonds = nullptr;
		plane_angles = nullptr;
		torsion = nullptr;
		moments = nullptr;
		rotations = nullptr;
		rotor = ASYMMETRIC;
		principle = nullptr;
		atoms = new std::vector<Atom>();
	}

	//Defined in molecule.cpp
	~DefaultMolecule();

	const std::vector<Atom> &getAtoms() const override {
		return (*(this->atoms));
	}

	void addAtom(Atom a) {
		this->atoms->push_back(a);
	}

	void setDistances(const Matrix<double> &dists) override {
		this->dists = new Matrix<double>(dists);
	}

	const Matrix<double> &getDistances() const override {
		return (*(this->dists));
	}

	void setBondAngles(const Matrix<double> &angles) override {
		this->bonds = new Matrix<double>(angles);
	}

	const Matrix<double> &getBondAngles() const override {
		return (*bonds);
	}

	void setPlaneAngles(const Matrix<double> &angles) override {
		this->plane_angles = new Matrix<double>(angles);
	}

	const Matrix<double> &getPlaneAngles() const override {
		return (*(this->plane_angles));
	}

	void setTorsionAngles(const Matrix<double> &angles) override {
		this->torsion = new Matrix<double>(angles);
	}

	const Matrix<double> &getTorsionAngles() const override {
		return (*(this->torsion));
	}

	void setMoments(const compchem::Matrix<double> &moms) override {
		this->moments = new compchem::Matrix<double>(moms);
	}

	const compchem::Matrix<double> &getMoments() const override {
		return (*(this->moments));
	}

	void setRotationalConstants(const std::vector<double> &rots) override {
		this->rotations = new std::vector<double>(rots);
	}

	const std::vector<double> &getRotationalConstants() const override {
		return (*(this->rotations));
	}

	void setRotorType(rotor_type rot) override {
		this->rotor = rot;
	}

	rotor_type getRotorType() const override {
		return (this->rotor);
	}

	void setPrincipleMoments(const std::vector<double> &moms) override {
		this->principle = new std::vector<double>(moms);
	}

	const std::vector<double> &getPrincipleMoments() const override {
		return (*(this->principle));
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
};

}
}



#endif /* MOLECULE_DEFAULT_HPP_ */
