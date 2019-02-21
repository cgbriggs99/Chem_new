/*
 * wavefunction_adapter.hpp
 *
 *  Created on: Feb 20, 2019
 *      Author: connor
 */

#ifndef WAVEFUNCTION_ADAPTER_HPP_
#define WAVEFUNCTION_ADAPTER_HPP_

#include "../Molecule/wavefunction.hpp"
#include "../Molecule/wavefunction_default.hpp"

namespace psi {
class Wavefunction;
class TwoElectronInt;
class MintsHelper;
}

namespace compchem {

class WavefunctionAdapter {
protected:
	WavefunctionAdapter() {
		;
	}
	static WavefunctionAdapter *singleton;
public:

	virtual ~WavefunctionAdapter() {
		;
	}

	static const WavefunctionAdapter *getSingelton() {
		if(singleton == nullptr) {
			singleton = new WavefunctionAdapter();
		}
		return (singleton);
	}

	static void init() {
		singleton = new WavefunctionAdapter();
	}

	//virtual compchem::AbstractWavefunction &convert(const psi::Wavefunction &wf, const psi::TwoElectronInt &tei) const;
	virtual compchem::AbstractWavefunction &convert(psi::MintsHelper &helper, double enuc) const;
};

}

#include "wavefunction_adapter.cpp"

#endif /* WAVEFUNCTION_ADAPTER_HPP_ */
