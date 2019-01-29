/*
 * math.cpp
 *
 *  Created on: Jan 28, 2019
 *      Author: cgbri
 */

#include "math.hpp"

double compchem::hypot3(double x, double y, double z) {
	double a, b, c;
	x = fabs(x);
	y = fabs(y);
	z = fabs(z);
	if(x >= y && x >= z) {
		a = x;
		b = y;
		c = z;
	} else if(y >= x && y >= z) {
		a = y;
		b = x;
		c = z;
	} else {
		a = z;
		b = x;
		c = y;
	}
	//Use this algorithm to reduce losses due to errors.
	return (a * sqrt(1 + (b * b) / (a * a) + (c * c) / (a * a)));
}

std::vector<double> compchem::crossprod(std::vector<double> v1, std::vector<double> v2) {
	std::vector<double> out;
	out[0] = v1[1] * v2[2] - v1[2] * v2[1];
	out[1] = v1[0] * v2[2] - v1[2] * v2[0];
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
	double off = diff * (a > b)? a: b;
	if(fabs(a - b) <= off) {
		return (0);
	} else if(a < b) {
		return (-1);
	} else {
		return (1);
	}
}
