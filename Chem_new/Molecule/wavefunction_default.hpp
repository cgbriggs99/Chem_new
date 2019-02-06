/*
 * wavefunction_default.hpp
 *
 *  Created on: Feb 4, 2019
 *      Author: Connor
 */

#ifndef WAVEFUNCTION_DEFAULT_HPP_
#define WAVEFUNCTION_DEFAULT_HPP_

#include "wavefunction.hpp"
#include "../Base/matrix_tei.hpp"
#include "molecule.hpp"

namespace compchem {
namespace strategies {

class DefaultWavefunction : public AbstractWavefunction {
private:
	compchem::Matrix<double> *_s, *_t, *_v, *_mux, *_muy, *_muz;
	double _enuc;
	compchem::strategies::TEIMatrix<double> *_tei;
	int n;
	int elecs;

public:
	DefaultWavefunction(const compchem::AbstractMolecule &mol) {
		_s = new compchem::Matrix<double>({mol.nelectron(), mol.nelectron()});
		_t = new compchem::Matrix<double>({mol.nelectron(), mol.nelectron()});
		_v = new compchem::Matrix<double>({mol.nelectron(), mol.nelectron()});
		_mux = new compchem::Matrix<double>({mol.nelectron(), mol.nelectron()});
		_muy = new compchem::Matrix<double>({mol.nelectron(), mol.nelectron()});
		_muz = new compchem::Matrix<double>({mol.nelectron(), mol.nelectron()});
		_enuc = 0;
		_tei = new compchem::strategies::TEIMatrix<double>(mol.nelectron());
		n = mol.norbital();
		elecs = mol.nelectron();
	}

	virtual ~DefaultWavefunction() {
		delete _s;
		delete _t;
		delete _v;
		delete _mux;
		delete _muy;
		delete _muz;
		delete _tei;
	}

	const compchem::AbstractMatrix<double> &s() const {
		return (*_s);
	}

	const compchem::AbstractMatrix<double> &t() const {
		return (*_t);
	}

	const compchem::AbstractMatrix<double> &v() const {
		return (*_v);
	}

	const compchem::AbstractMatrix<double> &two_electron() const {
		return (*_tei);
	}

	double enuc() const {
		return (_enuc);
	}

	const compchem::AbstractMatrix<double> &mux() const {
		return (*_mux);
	}

	const compchem::AbstractMatrix<double> &muy() const {
		return (*_muy);
	}

	const compchem::AbstractMatrix<double> &muz() const {
		return (*_muz);
	}

	void setS(double val, int i, int j) {
		_s->setEntry(val, i, j);
	}

	void setT(double val, int i, int j) {
		_t->setEntry(val, i, j);
	}

	void setV(double val, int i, int j) {
		_v->setEntry(val, i, j);
	}

	void setTEI(double val, int i, int j, int k, int l) {
		_tei->setEntry(val, i, j, k, l);
	}

	void setEnuc(double val) {
		_enuc = val;
	}

	void setMuX(double val, int i, int j) {
		_mux->setEntry(val, i, j);
	}

	void setMuY(double val, int i, int j) {
		_muy->setEntry(val, i, j);
	}

	void setMuZ(double val, int i, int j) {
		_muz->setEntry(val, i, j);
	}

	void setS(compchem::Matrix<double> *__restrict__ mat) {
		delete _s;
		_s = mat;
	}

	void setT(compchem::Matrix<double> *__restrict__ mat) {
		delete _t;
		_t = mat;
	}

	void setV(compchem::Matrix<double> *__restrict__ mat) {
		delete _v;
		_v = mat;
	}

	void setMuX(compchem::Matrix<double> *__restrict__ mat) {
		delete _mux;
		_mux = mat;
	}

	void setMuY(compchem::Matrix<double> *__restrict__ mat) {
		delete _muy;
		_muy = mat;
	}

	void setMuZ(compchem::Matrix<double> *__restrict__ mat) {
		delete _muz;
		_muz = mat;
	}

	void setTEI(compchem::strategies::TEIMatrix<double> *__restrict__ mat) {
		delete _tei;
		_tei = mat;
	}

	int getSize() const {
		return (this->n);
	}

	int nelectron() const {
		return (this->elecs);
	}
};


}
}



#endif /* WAVEFUNCTION_DEFAULT_HPP_ */
