/*
 * hessian_test.hpp
 *
 *  Created on: Feb 4, 2019
 *      Author: connor
 */

#include "test.hpp"
#include "../Project-2/hessian_default.hpp"
#include "../Base/matrix_default.hpp"
#include "../Base/eigenvalues_default.hpp"
#include "../Molecule/molecule_default.hpp"
#include "../Base/math.hpp"
#include <unistd.h>



class HessianTest : public test::Test {
private:
	compchem::EigenvalueStrategy<double> *strat;
	compchem::HessianStrategy *hess_strat;
	compchem::AbstractMolecule *mol;
	compchem::Matrix<double> *hessian;
	const char *dir;
public:
	HessianTest(const char *dir) {
		strat = new compchem::strategies::LapackEigenvalues<double>();
		hess_strat = new compchem::strategies::DefaultHessianStrategy();
		mol = new compchem::strategies::DefaultMolecule();
		hessian = nullptr;
		this->dir = dir;
	}

	~HessianTest() {
		delete strat;
		delete hess_strat;
		delete mol;
		if(hessian != nullptr) {
			delete hessian;
			hessian = nullptr;
		}
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

	void readHessianFile(const char *filename) {
		FILE *fp = fopen(filename, "r");
		int num = 0;
		int ignore = fscanf(fp, "%d", &num);
		assert(num == mol->natom());

		hessian = new compchem::Matrix<double>({3 * num, 3 * num});
		for(int i = 0; i < mol->natom() * 3; i++) {
			for(int j = 0; j < mol->natom(); j++) {
				double a, b, c;
				fscanf(fp, "%lf %lf %lf", &a, &b, &c);
				hessian->setEntry(a, i, 3 * j);
				hessian->setEntry(b, i, 3 * j + 1);
				hessian->setEntry(c, i, 3 * j + 2);
			}
		}
		fclose(fp);
	}

	void compareHessian(const compchem::AbstractMatrix<double> &eigs, const char *filename) {
		FILE *fp = fopen(filename, "r");
		while(!feof(fp)) {
			int n;
			double val;
			fscanf(fp, "%d %lf", &n, &val);
			assert(compchem::compareDoubles(eigs.getEntry(n), val, 0.01) == 0);
		}
	}

	void runTest() {
		chdir(dir);
		readGeomFile("geometry");
		readHessianFile("hessian");
		compchem::AbstractMatrix<double> &weight = hess_strat->massWeightMatrix(*mol, *hessian);
		compchem::Matrix<double> &eigs = strat->eigenvals(weight);

		eigs.sort();

		compchem::AbstractMatrix<double> &freqs = hess_strat->computeHessianFreqs(eigs);

		compareHessian(eigs, "eigenvalues");
		compareHessian(freqs, "frequencies");

		delete &weight;
		delete &eigs;
		delete &freqs;
		chdir("../");
	}
};

int main(void) {
	if(chdir("./data/hessian") == -1) {
		chdir("./Test/data/hessian");
	}

	HessianTest water("water"), benzene("benzene"), chlorobutene("3-chloro-1-butene");

	water.runTest();
	benzene.runTest();
	chlorobutene.runTest();

	return (0);
}
