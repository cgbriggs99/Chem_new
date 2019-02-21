/*
 * @BEGIN LICENSE
 *
 * plugin_test by Psi4 Developer, a plugin to:
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/liboptions/liboptions.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/matrix.h"
#include "../../Psi4-Bridging/wavefunction_adapter.hpp"
#include "../../Psi4-Bridging/matrix_adapter.hpp"
#include "../../Psi4-Bridging/tei_adapter.hpp"
#include "../../Project-3/scf_default.hpp"
#include "../../Base/base.hpp"
#include "../../Molecule/sto3g_basis_set.hpp"
#include "../../Base/eigenvalues_default_double.cpp"

namespace psi {
namespace plugin_test {

extern "C" PSI_API
int read_options(std::string name, Options& options) {
	if(name == "PLUGIN_TEST" || options.read_globals()) {
		/*- The amount of information printed to the output file -*/
		options.add_int("PRINT", 1);
	}

	return true;
}



extern "C" PSI_API
SharedWavefunction plugin_test(SharedWavefunction ref_wfn, Options& options) {
	int print = options.get_int("PRINT");

	/* Your code goes here */

	MintsHelper helper(ref_wfn);
	compchem::WavefunctionAdapter::init();
	compchem::TEIAdapter::init();
	compchem::Psi4MatrixAdapter::init();

	compchem::AbstractWavefunction *wfn =
	        &compchem::WavefunctionAdapter::getSingelton()->convert(helper,
	                ref_wfn.get()->molecule().get()->nuclear_repulsion_energy(
	                        std::array<double, 3>( {0, 0, 0})));

	double energy;
	compchem::strategies::DefaultSCFStrategy<
	        compchem::strategies::LapackEigenvalues<double>,
	        compchem::strategies::DefaultMatrixArithmeticStrategy<double> > scf;
	scf.runSCF(*wfn, nullptr, nullptr, nullptr, nullptr, &energy);

	ref_wfn.get()->set_efzc(energy);
	printf("%lf", energy);
	delete wfn;

	// Typically you would build a new wavefunction and populate it with data
	return ref_wfn;
}

}
} // End namespaces

