/*
 * scf_test.cpp
 *
 *  Created on: Feb 4, 2019
 *      Author: Connor
 */

#include "test.hpp"
#include "../Project-3/scf_default.hpp"
#include "../Base/matrix_tei.hpp"
#include "../Base/math.hpp"
#include "../Molecule/molecule_default.hpp"
#include "../Molecule/wavefunction_default.hpp"
#include "../Molecule/sto3g_basis_set.hpp"
#include "../Molecule/dz_basis_set.hpp"
#include <unistd.h>

template<typename _T>
class SCFTest : public test::Test {
private:
	compchem::SCFStrategy *strat;
	compchem::strategies::DefaultWavefunction<_T> *wfn;
	compchem::AbstractMolecule *mol;
	const char *dir;
public:
	SCFTest(const char *dir) {
		strat = new compchem::strategies::DefaultSCFStrategy();
		this->dir = dir;
		wfn = nullptr;
		mol = new compchem::strategies::DefaultMolecule();
	}

	~SCFTest() {
		delete strat;
		if(wfn != nullptr) {
			delete wfn;
		}
		delete mol;
	}

	void compare2d(const compchem::Matrix<double> &mat, const char *filename) {
		FILE *fp = fopen(filename, "r");

		for(int i = 0; i < mat.getShape(0); i++) {
			for(int j = 0; j < mat.getShape(1); j++) {
				double val = 0;
				int ignore = fscanf(fp, "%lf", &val);
				assert(compchem::compareDoubles(mat.getEntry(i, j), val, 0.001) == 0);
			}
		}
		fclose(fp);
	}

	void compareValue(double val, const char *filename) {
		FILE *fp = fopen(filename, "r");
		double expect = 0;
		int ignore = fscanf(fp, "%lf", &expect);
		assert(compchem::compareDoubles(val, expect, 0.001) == 0);
		fclose(fp);
	}

	void compareList(const std::vector<double> &vals, const char *filename) {
		FILE *fp = fopen(filename, "r");
		int i = 0;
		while(!feof(fp)) {
			double val = 0;
			int ignore = fscanf(fp, "%lf", &val);
			assert(compchem::compareDoubles(vals[i], val, 0.001) == 0);
			i++;
		}
		fclose(fp);
	}

	double readValueFile(const char *filename) {
		FILE *fp = fopen(filename, "r");
		double out = 0;
		int ignore = fscanf(fp, "%lf", &out);
		fclose(fp);
		return (out);
	}

	void readGeomFile(const char *filename) {
		int num = 0;
		FILE *fp = fopen(filename, "r");
		int ignore = fscanf(fp, "%d", &num);

		for (int i = 0; i < num; i++) {
			double n = 0, x = 0, y = 0, z = 0;
			ignore = fscanf(fp, "%lf %lf %lf %lf", &n, &x, &y, &z);
			mol->addAtom(compchem::Atom(n, compchem::amu(n), 0, x, y, z));
		}
		fclose(fp);
		return;
	}

	compchem::Matrix<double> &read2dFile(int n, const char *filename) {
		FILE *fp = fopen(filename, "r");

		compchem::Matrix<double> *out = new compchem::Matrix<double>({n, n});

		while(!feof(fp)) {
				double a = 0, b = 0, c = 0;
				fscanf(fp, "%lf %lf %lf", &a, &b, &c);
				out->setEntry(c, (int) a - 1, (int) b - 1);
		}
		fclose(fp);
		return (*out);
	}

	compchem::Matrix<double> &read2dSymmFile(int n, const char *filename) {
			FILE *fp = fopen(filename, "r");

			compchem::Matrix<double> *out = new compchem::Matrix<double>({n, n});

			while(!feof(fp)) {
					double a = 1, b = 1, c = 0;
					fscanf(fp, "%lf %lf %lf", &a, &b, &c);
					out->setEntry(c, (int) a - 1, (int) b - 1);
					out->setEntry(c, (int) b - 1, (int) a - 1);
			}
			fclose(fp);
			return (*out);
		}

	compchem::strategies::TEIMatrix<double> &read4dFile(int n, const char *filename) {
		FILE *fp = fopen(filename, "r");

		compchem::strategies::TEIMatrix<double> *out = new compchem::strategies::TEIMatrix<double>(n);

		while(!feof(fp)) {
			double a = 0, b = 0, c = 0, d = 0, e = 0;
			fscanf(fp, "%lf %lf %lf %lf %lf", &a, &b, &c, &d, &e);
			out->setEntry(e, (int) a - 1, (int) b - 1, (int) c - 1, (int) d - 1);
		}
		fclose(fp);
		return (*out);
	}

	void runTest() {
		chdir(dir);

		readGeomFile("geometry");
		//TODO Initialize the wavefunction.
		wfn = new compchem::strategies::DefaultWavefunction<_T>(*mol);
		wfn->setS(&read2dSymmFile(wfn->getSize(), "s"));
		wfn->setT(&read2dSymmFile(wfn->getSize(), "t"));
		wfn->setV(&read2dSymmFile(wfn->getSize(), "v"));
		wfn->setMuX(&read2dSymmFile(wfn->getSize(), "mux"));
		wfn->setMuY(&read2dSymmFile(wfn->getSize(), "muy"));
		wfn->setMuZ(&read2dSymmFile(wfn->getSize(), "muz"));

		wfn->setEnuc(readValueFile("enuc"));
		wfn->setTEI(&read4dFile(wfn->getSize(), "eri"));

		compchem::Matrix<double> *hamiltonian = (compchem::Matrix<double> *) &strat->findHamiltonian(*wfn);
		compchem::Matrix<double> *fock, *c, *density;
		double energy;
		strat->runSCF(*wfn, (compchem::AbstractMatrix<double> **) &fock, (compchem::AbstractMatrix<double> **) &c,
				(compchem::AbstractMatrix<double> **) &density, &energy);
		std::vector<double> *charges = &strat->findElectronCharge(*mol, *wfn, *density);
		std::array<double, 3> *moment = &strat->findDipole(*mol, *density, *wfn);

		compare2d(*density, "density");
		compareList(*charges, "charges");
		compareValue(energy, "etotal");
		compare2d(*hamiltonian, "hamiltonian");
		compareList(std::vector<double>(moment->begin(), moment->end()), "moment");

		delete hamiltonian;
		delete fock;
		delete c;
		delete density;
		delete charges;
		delete moment;
		chdir("../");
	}

};


int main(void) {
	if(chdir("./data/scf") == -1) {
		chdir("./Test/data/scf");
	}

	SCFTest<compchem::strategies::STO3GBasisSet> sto3g_water("sto3g-water"), sto3g_methane("sto3g-methane");
	SCFTest<compchem::strategies::DZBasisSet> dz_water("dz-water");

	sto3g_water.runTest();
	dz_water.runTest();
	sto3g_methane.runTest();
	return (0);
}

