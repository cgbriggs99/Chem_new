#include "test.hpp"
#include "../Project-3/scf_default.hpp"
#include "../Project-8/diis_scf.hpp"
#include "../Project-12/cis_default.hpp"
#include "../Molecule/spin_adapt_default.hpp"
#include "../Base/base.hpp"
#include "../Molecule/molecule_default.hpp"
#include "../Molecule/wavefunction_default.hpp"
#include "../Molecule/sto3g_basis_set.hpp"
#include "../Molecule/dz_basis_set.hpp"
#include <unistd.h>


template<typename _T>
class CISTest : public test::Test {
private:
	compchem::AbstractCIEnergyStrategy *strat;
	compchem::SCFStrategy *scf;
	compchem::strategies::DefaultWavefunction<_T> *wfn;
	compchem::AbstractMolecule *mol;
	compchem::AbstractSpinAdaptationStrategy *spin;
	const char *dir;
public:
	CISTest(const char *dir) {
		strat = new compchem::strategies::DefaultCISEnergyStrategy<
		        compchem::strategies::LapackEigenvalues<double> >();
		scf =
		        new compchem::strategies::SCF_DIISStrategy<
		                compchem::strategies::LapackEigenvalues<double>,
		                compchem::strategies::DefaultMatrixArithmeticStrategy<
		                        double> >();
		this->dir = dir;
		wfn = nullptr;
		mol = new compchem::strategies::DefaultMolecule();
		spin = new compchem::strategies::DefaultSpinAdaptationStrategy();
	}

	~CISTest() {
		delete strat;
		if(wfn != nullptr) {
			delete wfn;
		}
		delete mol;
		delete scf;
		delete spin;
	}

	void compare2d(const compchem::Matrix<double> &mat, const char *filename) {
		FILE *fp = fopen(filename, "r");

		for(int i = 0; i < mat.getShape(0); i++) {
			for(int j = 0; j < mat.getShape(1); j++) {
				double val = 0;
				int ignore = fscanf(fp, "%lf", &val);
				assert(
				        compchem::compareDoubles(mat.getEntry(i, j), val, 0.001)
				                == 0);
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
			if(ignore <= 0) {
				break;
			}
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

		for(int i = 0; i < num; i++) {
			double n = 0, x = 0, y = 0, z = 0;
			ignore = fscanf(fp, "%lf %lf %lf %lf", &n, &x, &y, &z);
			mol->addAtom(compchem::Atom(n, compchem::amu(n), 0, x, y, z));
		}
		fclose(fp);
		return;
	}

	compchem::Matrix<double> &read2dFile(int n, const char *filename) {
		FILE *fp = fopen(filename, "r");

		compchem::Matrix<double> *out = new compchem::Matrix<double>( {n, n});

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

		compchem::Matrix<double> *out = new compchem::Matrix<double>( {n, n});

		while(!feof(fp)) {
			double a = 1, b = 1, c = 0;
			fscanf(fp, "%lf %lf %lf", &a, &b, &c);
			out->setEntry(c, (int) a - 1, (int) b - 1);
			out->setEntry(c, (int) b - 1, (int) a - 1);
		}
		fclose(fp);
		return (*out);
	}

	compchem::strategies::TEIMatrix<double> &read4dFile(int n,
	        const char *filename) {
		FILE *fp = fopen(filename, "r");

		compchem::strategies::TEIMatrix<double> *out =
		        new compchem::strategies::TEIMatrix<double>(n);

		int lines = 0;

		while(!feof(fp)) {
			double a = 0, b = 0, c = 0, d = 0, e = 0;
			int count = fscanf(fp, "%lf %lf %lf %lf %lf", &a, &b, &c, &d, &e);
			lines++;
			if(count <= 0) {
				continue;
			}

			out->setEntry(e, (int) a - 1, (int) b - 1, (int) c - 1,
			        (int) d - 1);
		}
		fclose(fp);
		return (*out);
	}

	void runTest() {
		chdir(dir);

		readGeomFile("geometry");
		wfn = new compchem::strategies::DefaultWavefunction<_T>(*mol);
		wfn->setS(&read2dSymmFile(wfn->getSize(), "s"));
		wfn->setT(&read2dSymmFile(wfn->getSize(), "t"));
		wfn->setV(&read2dSymmFile(wfn->getSize(), "v"));
		wfn->setEnuc(readValueFile("enuc"));
		wfn->setTEI(&read4dFile(wfn->getSize(), "eri"));

		compchem::Matrix<double> *hamiltonian =
		        (compchem::Matrix<double> *) &scf->findHamiltonian(*wfn), *fock,
		        *c, *eigs;
		double energy;
		scf->runSCF(*wfn, (compchem::AbstractMatrix<double> **) &fock,
		        (compchem::AbstractMatrix<double> **) &c, nullptr,
		        (compchem::AbstractMatrix<double> **) &eigs, &energy);

		delete hamiltonian;
		hamiltonian = (compchem::Matrix<double> *) &scf->findHamiltonian(*wfn);

		compchem::Matrix<double> *sofock = (compchem::Matrix<double> *) &spin
		        ->spinAdapt1elec(*fock, *c);
		compchem::Matrix<double> *sotei = (compchem::Matrix<double> *) &spin
		        ->spinAdaptTEI(wfn->two_electron(), *c);

		//Make sure stuff is done properly.
		double sum = 0;
		for(int i = 0; i < mol->nelectron(); i++) {
			for(int j = 0; j < mol->nelectron(); j++) {
				for(int a = mol->nelectron(); a < sotei->getShape(0); a++) {
					for(int b = mol->nelectron(); b < sotei->getShape(0); b++) {
						sum += (sotei->getEntry(i, j, a, b)
						        * sotei->getEntry(i, j, a, b))
						        / (sofock->getEntry(i, i)
						                + sofock->getEntry(j, j)
						                - sofock->getEntry(a, a)
						                - sofock->getEntry(b, b));
					}
				}
			}
		}

		compareValue(sum / 4, "mp2_energy");

		compchem::Matrix<double> *cis_energy =
		        (compchem::Matrix<double> *) &strat->CIEnergy(*sofock, *sotei,
		                mol->nelectron());

		std::vector<double> vec;
		for(int i = 0; i < cis_energy->getShape(0); i++) {
			vec.push_back(cis_energy->getEntry(i));
		}

		compareList(vec, "cis_energy");

		delete eigs;
		delete fock;
		delete c;
		delete hamiltonian;
		delete cis_energy;
		delete sofock;
		delete sotei;
		chdir("../");
	}

};

int main(void) {
	if(chdir("./data/energies") == -1) {
		chdir("./Test/data/energies");
	}

	CISTest<compchem::strategies::STO3GBasisSet> sto3g_water("sto3g-water");
	sto3g_water.runTest();
	CISTest<compchem::strategies::STO3GBasisSet> sto3g_methane("sto3g-methane");
	sto3g_methane.runTest();
	CISTest<compchem::strategies::DZBasisSet> dz_water("dz-water");
	dz_water.runTest();

	return (0);
}
