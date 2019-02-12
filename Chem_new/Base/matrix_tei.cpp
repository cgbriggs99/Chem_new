/*
 * matrix_tei.cpp
 *
 *  Created on: Feb 4, 2019
 *      Author: Connor
 */

#ifndef __MATRIX_TEI_CPP__
#define __MATRIX_TEI_CPP__

#include "matrix_tei.hpp"
#include <stdarg.h>
#include <stdlib.h>

static int tei_index(int a, int b, int c, int d) {
	int ab, cd;
	if (a > b) {
		ab = (a * (a + 1)) / 2 + b;
	} else {
		ab = (b * (b + 1)) / 2 + a;
	}

	if (c > d) {
		cd = (c * (c + 1)) / 2 + d;
	} else {
		cd = (d * (d + 1)) / 2 + c;
	}

	if (ab > cd) {
		return ((ab * (ab + 1)) / 2 + cd);
	} else {
		return ((cd * (cd + 1)) / 2 + ab);
	}
}

template<typename _T, typename _Alloc>
compchem::strategies::TEIMatrix<_T, _Alloc>::TEIMatrix(int n) {
	this->data = new _T[tei_index(n, n, n, n)];
	this->n = n;
	this->size = n * n * n * n;
}

template<typename _T, typename _Alloc>
compchem::strategies::TEIMatrix<_T, _Alloc>::~TEIMatrix() {
	delete[] data;
}

template<typename _T, typename _Alloc>
const _T &compchem::strategies::TEIMatrix<_T, _Alloc>::getEntry(int index,
		...) const {
	va_list list;
	va_start(list, index);

	int a, b, c, d;
	a = index;
	b = va_arg(list, int);
	c = va_arg(list, int);
	d = va_arg(list, int);
	va_end(list);
	return (this->data[tei_index(a, b, c, d)]);
}

template<typename _T, typename _Alloc>
const _T &compchem::strategies::TEIMatrix<_T, _Alloc>::getEntry(
		std::vector<int> index) const {
	return (this->data[tei_index(index[0], index[1], index[2], index[3])]);
}

template<typename _T, typename _Alloc>
void compchem::strategies::TEIMatrix<_T, _Alloc>::setEntry(const _T &ent, int index,
		...) {
	va_list list;
	va_start(list, index);

	int a, b, c, d;
	a = index;
	b = va_arg(list, int);
	c = va_arg(list, int);
	d = va_arg(list, int);
	va_end(list);
	this->data[tei_index(a, b, c, d)] = ent;
}

template<typename _T, typename _Alloc>
void compchem::strategies::TEIMatrix<_T, _Alloc>::setEntry(const _T &ent, std::vector<int> index) {
	this->data[tei_index(index[0], index[1], index[2], index[3])] = ent;
}

template<typename _T, typename _Alloc>
void compchem::strategies::TEIMatrix<_T, _Alloc>::setEntry(const _T &&ent, int index,
		...) {
	va_list list;
	va_start(list, index);

	int a, b, c, d;
	a = index;
	b = va_arg(list, int);
	c = va_arg(list, int);
	d = va_arg(list, int);
	va_end(list);
	this->data[tei_index(a, b, c, d)] = ent;
}

template<typename _T, typename _Alloc>
void compchem::strategies::TEIMatrix<_T, _Alloc>::setEntry(const _T &&ent, std::vector<int> index) {
	this->data[tei_index(index[0], index[1], index[2], index[3])] = ent;
}

template<typename _T, typename _Alloc>
int compchem::strategies::TEIMatrix<_T, _Alloc>::getSize() const {
	return (this->size);
}

template<typename _T, typename _Alloc>
int compchem::strategies::TEIMatrix<_T, _Alloc>::getShape(int dim) const {
	if(dim < 4 && dim >= 0) {
		return (this->n);
	} else {
		return (0);
	}
}

template<typename _T, typename _Alloc>
int compchem::strategies::TEIMatrix<_T, _Alloc>::getDimension() const {
	return (4);
}

#endif
