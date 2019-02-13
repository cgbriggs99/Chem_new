/*
 * math.cpp
 *
 *  Created on: Jan 28, 2019
 *      Author: cgbri
 */

#include "math.hpp"
#include <cmath>
#include <exception>
#include <stdexcept>

double compchem::hypot3(double x, double y, double z) {
	double a, b, c, x1, y1, z1;
	x1 = fabs(x);
	y1 = fabs(y);
	z1 = fabs(z);
	if(x1 >= y1 && x1 >= z1) {
		a = x1;
		b = y1;
		c = z1;
	} else if(y1 >= x1 && y1 >= z1) {
		a = y1;
		b = x1;
		c = z1;
	} else {
		a = z1;
		b = x1;
		c = y1;
	}
	//Use this algorithm to reduce losses due to errors.
	return (a * sqrt(1 + (b * b) / (a * a) + (c * c) / (a * a)));
}

std::vector<double> compchem::crossprod(std::vector<double> v1, std::vector<double> v2) {
	std::vector<double> out(3);
	out[0] = v1[1] * v2[2] - v1[2] * v2[1];
	out[1] = v1[2] * v2[0] - v1[0] * v2[2];
	out[2] = v1[0] * v2[1] - v1[1] * v2[0];
	return (out);
}

double compchem::dotprod(std::vector<double> v1, std::vector<double> v2) {
	double out = 0;
	for(int i = 0; i < v1.size(); i++) {
		out += v1[i] * v2[i];
	}
	return (out);
}

int compchem::compareDoubles(double a, double b, double diff) {
	double off = diff * ((fabs(a) > fabs(b))? fabs(a): fabs(b));
	if(std::isnan(fabs(a)) || std::isnan(fabs(b)) || std::isinf(fabs(a)) || std::isinf(fabs(b))) {
		throw(new std::domain_error("Testing infinity or NaN!"));
	}
	if(fabs(a - b) <= diff || fabs(a - b) <= off) {
		return (0);
	} else if(a < b) {
		return (-1);
	} else {
		return (1);
	}
}
