/*
 * basis_set.hpp
 *
 *  Created on: Feb 4, 2019
 *      Author: Connor
 */

#ifndef WAVEFUNCTION_HPP_
#define WAVEFUNCTION_HPP_

#include "../Base/matrix.hpp"
#include "../Base/matrix_default.hpp"

namespace compchem {

class AbstractWavefunction {
public:
	AbstractWavefunction() {
		;
	}
	virtual ~AbstractWavefunction() = default;

	virtual const compchem::AbstractMatrix<double> &s() const = 0;
	virtual const compchem::AbstractMatrix<double> &t() const = 0;
	virtual const compchem::AbstractMatrix<double> &v() const = 0;
	virtual const compchem::AbstractMatrix<double> &two_electron() const = 0;
	virtual double enuc() const = 0;
	virtual const compchem::AbstractMatrix<double> &mux() const = 0;
	virtual const compchem::AbstractMatrix<double> &muy() const = 0;
	virtual const compchem::AbstractMatrix<double> &muz() const = 0;
	virtual int getSize() const = 0;
	virtual int nelectron() const = 0;
};

}

#endif /* WAVEFUNCTION_HPP_ */
