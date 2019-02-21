/*
 * tei_adapter.hpp
 *
 *  Created on: Feb 20, 2019
 *      Author: connor
 */

#ifndef TEI_ADAPTER_HPP_
#define TEI_ADAPTER_HPP_

#include "../Base/base.hpp"

namespace psi {
class TwoElectronInt;
class Matrix;
typedef std::shared_ptr<Matrix> SharedMatrix;
}

namespace compchem {

class TEIAdapter {
protected:
	TEIAdapter() {
		;
	}

	static TEIAdapter *singleton;
public:

	virtual ~TEIAdapter() {
		;
	}

	static const TEIAdapter *getSingleton() {
		if(singleton == nullptr) {
			singleton = new TEIAdapter();
		}
		return (singleton);
	}

	static void init() {
		singleton = new TEIAdapter();
	}

	virtual compchem::strategies::TEIMatrix<double> &convert(psi::SharedMatrix tei) const;

};

}

#include "tei_adapter.cpp"

#endif /* TEI_ADAPTER_HPP_ */
