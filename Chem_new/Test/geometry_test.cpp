/*
 * geometry_test.cpp
 *
 *  Created on: Jan 28, 2019
 *      Author: connor
 */

#include "../Project-1/geometry.hpp"
#include "/usr/local/psi4/include/psi4/pragma.h"

#include "test.hpp"
#include "/usr/local/psi4/include/psi4/libmints/molecule.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "../Project-1/geometry_default.hpp"
#include "../Base/math.hpp"
#include <string.h>

class GeometryTest : public test::Test {
private:
	psi::Molecule *molecule;
	compchem::GeometryCalcStrategy *strat;
	compchem::Matrix<double> *dists, *bonds, *planes, *torsions, *moment;
	std::vector<double> *principles;
	std::vector<double> *com;
	const char *dir;
public:
	GeometryTest(const char *dir) {
		molecule = new psi::Molecule();
		strat = new compchem::strategies::DefaultGeometryStrategy();
		this->dir = dir;
		dists = nullptr;
		bonds = nullptr;
		planes = nullptr;
		torsions = nullptr;
		moment = nullptr;
		principles = nullptr;
		com = nullptr;
	}

	~GeometryTest() {
		delete molecule;
		delete strat;
		delete dists;
		delete bonds;
		delete planes;
		delete torsions;
		delete moment;
		delete principles;
		delete com;
	}

	void readGeomFile(const char *filename) {
		int num;
		FILE *fp = fopen(filename, "r+");
		int ignore = fscanf(fp, "%d", &num);



		for(int i = 0; i < num; i++) {
			double n, x, y, z;
			ignore = fscanf(fp, "%lf %lf %lf %lf", &n, &x, &y, &z);
			molecule->add_atom(n, x, y, z);
		}
		fclose(fp);
		return;
	}

	void compareDistances(const char *filename) {
		FILE *fp = fopen(filename, "r");
		char buffer[1000];

		while(!feof(fp)) {
			fgets(buffer, 998, fp);
			int a, b;
			double dist;
			sscanf(buffer, "%d %d %lf", &a, &b, &dist);
			assert(compchem::compareDoubles(dists->getEntry({a, b}), dist, 0.0001) == 0);
		}
		fclose(fp);
	}

	void compareBondAngles(const char *filename) {
		FILE *fp = fopen(filename, "r");
		char buffer[1000];

		while(!feof(fp)) {
			fgets(buffer, 998, fp);
			int a, b, c;
			double dist;
			sscanf(buffer, "%d- %d- %d %lf", &a, &b, &c, &dist);
			assert(compchem::compareDoubles(bonds->getEntry({a, b, c}) * 180 / M_PI, dist, 0.01) == 0);
		}
		fclose(fp);
	}

	void comparePlaneAngles(const char *filename) {
		FILE *fp = fopen(filename, "r");
		char buffer[1000];

		while(!feof(fp)) {
			fgets(buffer, 998, fp);
			int a, b, c, d;
			double dist;
			int scans = sscanf(buffer, "%d- %d- %d- %d %lf", &a, &b, &c, &d, &dist);
			assert(compchem::compareDoubles(planes->getEntry({a, b, c, d}) * 180 / M_PI, dist, 0.0001) == 0);
		}
		fclose(fp);
	}

	void compareTorsionAngles(const char *filename) {
		FILE *fp = fopen(filename, "r");
		char buffer[1000];

		while(!feof(fp)) {
			fgets(buffer, 998, fp);
			int a, b, c, d;
			double dist;
			sscanf(buffer, "%d- %d- %d- %d %lf", &a, &b, &c, &d, &dist);
			assert(compchem::compareDoubles(torsions->getEntry({a, b, c, d}) * 180 / M_PI, dist, 0.0001) == 0);
		}
		fclose(fp);
	}

	void compareCOM(const char *filename, const std::vector<double> &com) {
		FILE *fp = fopen(filename, "r");
		double x, y, z;
		fscanf(fp, "%lf %lf %lf", &x, &y, &z);
		assert(compchem::compareDoubles(com[0], x, 0.01) == 0);
		assert(compchem::compareDoubles(com[1], y, 0.01) == 0);
		assert(compchem::compareDoubles(com[2], z, 0.01) == 0);
		fclose(fp);
	}

	void compareMoments(const char *filename) {
		FILE *fp = fopen(filename, "r");
		double x, y, z;
		for(int i = 0; i < 3; i++) {
			fscanf(fp, "%lf %lf %lf", &x, &y, &z);
			assert(compchem::compareDoubles(moment->getEntry({i, 0}), x, 0.001) == 0);
			assert(compchem::compareDoubles(moment->getEntry({i, 1}), y, 0.001) == 0);
			assert(compchem::compareDoubles(moment->getEntry({i, 2}), z, 0.001) == 0);
		}
		fclose(fp);
	}

	void comparePrincipleMoments(const char *filename) {
		FILE *fp = fopen(filename, "r");
		double x, y, z;
		fscanf(fp, "%lf %lf %lf", &x, &y, &z);
		double a, b, c;

		if(x >= y && x >= z) {
			a = x;
			if(y > z) {
				b = y;
				c = z;
			} else {
				b = z;
				c = y;
			}
		} else if(y >= x && y >= z) {
			a = y;
			if(x > z) {
				b = x;
				c = z;
			} else {
				b = z;
				c = x;
			}
		} else {
			a = z;
			if(x > y) {
				b = x;
				c = y;
			} else {
				b = y;
				c = x;
			}
		}

		assert(compchem::compareDoubles(principles->at(0), c, 0.001) == 0);
		assert(compchem::compareDoubles(principles->at(1), b, 0.001) == 0);
		assert(compchem::compareDoubles(principles->at(2), a, 0.001) == 0);
		fclose(fp);
	}

	void runTest() {
		chdir(dir);
		readGeomFile("input");
		*dists = strat->findDistances(*molecule);
		*bonds = strat->findBondAngles(*molecule);
		*planes = strat->findPlaneAngles(*molecule, *bonds);
		*torsions = strat->findTorsionAngles(*molecule, *bonds);
		*com = strat->findCenterOfMass(*molecule);
		molecule->move_to_com();
		*moment = strat->findMoments(*molecule);
		*principles = strat->findPrincipleMoments(*moment);

		compareDistances("./distances");
		compareBondAngles("./bond_angles");
		compareCOM("./center_of_mass", *com);
		compareMoments("./moment_of_inertia");
		comparePlaneAngles("./plane_angles");
		compareTorsionAngles("./torsion_angles");
		comparePrincipleMoments("./principle_moments");

		chdir("..");
	}
};

int main(void) {
	chdir("./Test/data/geometry");
	GeometryTest acetaldehyde("./acetaldehyde"), benzene("./benzene"), allene("./allene");

	acetaldehyde.runTest();
	benzene.runTest();
	allene.runTest();
	return (0);
}
