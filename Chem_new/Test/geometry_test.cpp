/*
 * geometry_test.cpp
 *
 *  Created on: Jan 28, 2019
 *      Author: connor
 */

#include "../Project-1/geometry.hpp"

#include "test.hpp"
#include "../Molecule/molecule_default.hpp"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "../Project-1/geometry_default.hpp"
#include "../Base/math.hpp"
#include <string.h>

class GeometryTest : public test::Test {
private:
	compchem::AbstractMolecule *molecule;
	compchem::GeometryCalcStrategy *strat;
	const char *dir;
public:
	GeometryTest(const char *dir) {
		molecule = new compchem::strategies::DefaultMolecule();
		strat = new compchem::strategies::DefaultGeometryStrategy();
		this->dir = dir;
	}

	~GeometryTest() {
		delete molecule;
		delete strat;
	}

	void readGeomFile(const char *filename) {
		int num = 0;
		FILE *fp = fopen(filename, "r");
		int ignore = fscanf(fp, "%d", &num);



		for(int i = 0; i < num; i++) {
			double n, x, y, z;
			ignore = fscanf(fp, "%lf %lf %lf %lf", &n, &x, &y, &z);
			molecule->addAtom(compchem::Atom(n, compchem::amu(n), 0, x, y, z));
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
			assert(compchem::compareDoubles(molecule->getDistances().getEntry({a, b}), dist, 0.0001) == 0);
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
			assert(compchem::compareDoubles(molecule->getBondAngles().getEntry({a, b, c}) * 180 / M_PI, dist, 0.01) == 0);
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
			assert(compchem::compareDoubles(molecule->getPlaneAngles().getEntry({a, b, c, d}) * 180 / M_PI, dist, 0.0001) == 0);
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
			assert(compchem::compareDoubles(molecule->getTorsionAngles().getEntry({a, b, c, d}) * 180 / M_PI, dist, 0.0001) == 0);
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
			assert(compchem::compareDoubles(molecule->getMoments().getEntry({i, 0}), x, 0.001) == 0);
			assert(compchem::compareDoubles(molecule->getMoments().getEntry({i, 1}), y, 0.001) == 0);
			assert(compchem::compareDoubles(molecule->getMoments().getEntry({i, 2}), z, 0.001) == 0);
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

		assert(compchem::compareDoubles(molecule->getPrincipleMoments()[0], c, 0.001) == 0);
		assert(compchem::compareDoubles(molecule->getPrincipleMoments()[1], b, 0.001) == 0);
		assert(compchem::compareDoubles(molecule->getPrincipleMoments()[2], a, 0.001) == 0);
		fclose(fp);
	}

	void runTest() {
		chdir(dir);
		readGeomFile("input");
		compchem::Matrix<double> &dist = strat->findDistances(*molecule);
		molecule->setDistances(dist);
		compchem::Matrix<double> &bonds = strat->findBondAngles(*molecule);
		molecule->setBondAngles(bonds);
		compchem::Matrix<double> &plane = strat->findPlaneAngles(*molecule);
		compchem::Matrix<double> &torque = strat->findTorsionAngles(*molecule);
		molecule->setPlaneAngles(plane);
		molecule->setTorsionAngles(torque);
		std::vector<double> &com = strat->findCenterOfMass(*molecule);
		molecule->translateCOM(com);
		compchem::Matrix<double> &moms = strat->findMoments(*molecule);
		molecule->setMoments(moms);
		std::vector<double> &rots = strat->findPrincipleMoments(*molecule);
		molecule->setPrincipleMoments(rots);

		compareDistances("./distances");
		compareBondAngles("./bond_angles");
		compareCOM("./center_of_mass", com);
		compareMoments("./moment_of_inertia");
		comparePlaneAngles("./plane_angles");
		compareTorsionAngles("./torsion_angles");
		comparePrincipleMoments("./principle_moments");


		delete &dist;
		delete &bonds;
		delete &plane;
		delete &torque;
		delete &com;
		delete &moms;
		delete &rots;
		chdir("..");
	}
};

int main(void) {
	char buff[1000];
	 getcwd(buff, 999);
	 puts(buff);

	chdir("./data/geometry");
	GeometryTest acetaldehyde("./acetaldehyde"), benzene("./benzene"), allene("./allene");

	acetaldehyde.runTest();
//	benzene.runTest();
//	allene.runTest();
	return (0);
}
