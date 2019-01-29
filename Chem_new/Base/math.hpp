/*
 * math.hpp
 *
 *  Created on: Jan 28, 2019
 *      Author: cgbri
 */

#ifndef BASE_MATH_HPP_
#define BASE_MATH_HPP_

#include "matrix.hpp"

#include <math.h>
#include <vector>

namespace compchem {


/*
 * Functions to calculate these values for small vectors.
 */
double hypot3(double x, double y, double z);

std::vector<double> crossprod(std::vector<double> v1, std::vector<double> v2);
double dotprod(std::vector<double> v1, std::vector<double> v2);

//Compares a and b to within fraction_difference * (greater of a and b)
int compareDoubles(double a, double b, double fraction_difference);

}

#endif /* BASE_MATH_HPP_ */
