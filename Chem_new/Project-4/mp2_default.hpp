/*
 * mp2_default.hpp
 *
 *  Created on: Feb 11, 2019
 *      Author: connor
 */

#ifndef MP2_DEFAULT_HPP_
#define MP2_DEFAULT_HPP_


#include "mp2.hpp"

namespace compchem {
namespace strategies {

template<typename _T>
class DefaultMP2Strategy : public AbstractMP2Strategy<_T> {
private:
	_T *basis_set;
public:
	DefaultMP2Strategy() {
		basis_set = new _T();
	}
	virtual ~DefaultMP2Strategy() {
		delete basis_set;
	}

	double mp2Energy(const compchem::AbstractMolecule &mol, const compchem::AbstractMatrix<double> &tei,
			const compchem::AbstractMatrix<double> &orbs, const compchem::AbstractMatrix<double> &eigs) override;

};


}
}

#include "mp2_default.cpp"


#endif /* MP2_DEFAULT_HPP_ */
