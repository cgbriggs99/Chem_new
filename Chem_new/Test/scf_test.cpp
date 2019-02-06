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
#include <unistd.h>

class SCFTest : public test::Test {
private:
	compchem::SCFStrategy *strat;
	compchem::strategies::DefaultWavefunction *wfn;
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
				double val;
				int ignore = fscanf(fp, "%lf", &val);
				assert(compchem::compareDoubles(mat.getEntry(i, j), val, 0.001) == 0);
			}
		}
		fclose(fp);
	}

	void compareValue(double val, const char *filename) {
		FILE *fp = fopen(filename, "r");
		double expect;
		int ignore = fscanf(fp, "%lf", &expect);
		assert(compchem::compareDoubles(val, expect, 0.001) == 0);
		fclose(fp);
	}

	void compareList(const std::vector<double> &vals, const char *filename) {
		FILE *fp = fopen(filename, "r");
		int i = 0;
		while(!feof(fp)) {
			double val;
			int ignore = fscanf(fp, "%lf", &val);
			assert(compchem::compareDoubles(vals[i], val, 0.001) == 0);
			i++;
		}
		fclose(fp);
	}

	double readValueFile(const char *filename) {
		FILE *fp = fopen(filename, "r");
		double out;
		int ignore = fscanf(fp, "%lf", &out);
		fclose(fp);
		return (out);
	}

	void readGeomFile(const char *filename) {
		int num = 0;
		FILE *fp = fopen(filename, "r");
		int ignore = fscanf(fp, "%d", &num);

		for (int i = 0; i < num; i++) {
			double n, x, y, z;
			ignore = fscanf(fp, "%lf %lf %lf %lf", &n, &x, &y, &z);
			mol->addAtom(compchem::Atom(n, compchem::amu(n), 0, x, y, z));
		}
		fclose(fp);
		return;
	}

	compchem::Matrix<double> &read2dFile(int n, const char *filename) {
		FILE *fp = fopen(filename, "r");
		int num = 0;
		int ignore = fscanf(fp, "%d", &num);
		assert(num == mol->natom());

		compchem::Matrix<double> *out = new compchem::Matrix<double>({n, n});

		while(!feof(fp)) {
				double a, b, c;
				fscanf(fp, "%lf %lf %lf", a, b, c);
				out->setEntry(c, (int) a, (int) b);
		}
		fclose(fp);
		return (*out);
	}

	compchem::strategies::TEIMatrix<double> &read4dFile(int n, const char *filename) {
		FILE *fp = fopen(filename, "r");
		int num = 0;
		int ignore = fscanf(fp, "%d", &num);
		assert(num == mol->natom());

		compchem::strategies::TEIMatrix<double> *out = new compchem::strategies::TEIMatrix<double>(n);

		while(!feof(fp)) {
			double a, b, c, d, e;
			fscanf(fp, "%lf %lf %lf %lf %lf", a, b, c, d, e);
			out->setEntry(e, (int) a, (int) b, (int) c, (int) d);
		}
		fclose(fp);
		return (*out);
	}

	void runTest() {
		chdir(dir);

		readGeomFile("geometry");
		//TODO Initialize the wavefunction.
		wfn = new compchem::strategies::DefaultWavefunction(*mol);
		wfn->setS(&read2dFile(mol->norbital(), "s"));
		wfn->setT(&read2dFile(mol->norbital(), "t"));
		wfn->setV(&read2dFile(mol->norbital(), "v"));
		wfn->setMuX(&read2dFile(mol->norbital(), "mux"));
		wfn->setMuY(&read2dFile(mol->norbital(), "muy"));
		wfn->setMuZ(&read2dFile(mol->norbital(), "muz"));

		double enuc = readValueFile("enuc");
		wfn->setTEI(&read4dFile(mol->norbital(), "eri"));

		compchem::Matrix<double> *hamiltonian = (compchem::Matrix<double> *) &strat->findHamiltonian(*wfn);
		compchem::Matrix<double> *fock, *c, *density;
		double energy;
		strat->runSCF(*wfn, (compchem::AbstractMatrix<double> **) &fock, (compchem::AbstractMatrix<double> **) &c,
				(compchem::AbstractMatrix<double> **) &density, &energy);
		std::vector<double> *charges = &strat->findElectronCharge(*mol, *wfn, *density);
		std::array<double, 3> *moment = &strat->findDipole(*density, *wfn);

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

	SCFTest sto3g_water("sto3g-water");

	sto3g_water.runTest();
	return (0);
}

