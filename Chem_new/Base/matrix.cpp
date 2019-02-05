#ifndef __MATRIX_CPP__
#define __MATRIX_CPP__

#include "matrix_default.hpp"
#include <exception>
#include <stdexcept>
#include <string>
#include <initializer_list>
#include <stdlib.h>

template<typename _T, typename _Alloc>
compchem::Matrix<_T, _Alloc>::Matrix(std::initializer_list<int> shape) {
	this->dimensions = shape.size();
	this->shape = new int[this->dimensions];
	const int *list = shape.begin();
	for (int i = 0; i < this->dimensions; i++) {
		this->shape[i] = list[i];
	}
	int prod = 1;
	for (int i = 0; i < dimensions; i++) {
		prod *= this->shape[i];
	}
	this->data = this->allocator.allocate(prod);
	this->size = prod;
}

template<typename _T, typename _Alloc>
compchem::Matrix<_T, _Alloc>::Matrix(const AbstractMatrix<_T> &copy) {
	printf("Copied!\n");
	this->dimensions = copy.getDimension();
	this->size = copy.getSize();
	this->shape = new int[this->dimensions];
	for (int i = 0; i < this->dimensions; i++) {
		this->shape[i] = copy.getShape(i);
	}
	this->data = this->allocator.allocate(this->size);
	std::vector<int> index(this->dimensions);
	for (int i = 0; i < this->size; i++) {
		int hold = i;
		//TODO optimize this.
		for(int j = this->dimensions - 1; j >= 0; j--) {
			index[j] = hold % this->getShape(j);
			hold /= this->getShape(j);
		}
		this->data[i] = copy.getEntry(index);
	}
}

template<typename _T, typename _Alloc>
template<typename _U, typename _Alloc2>
compchem::Matrix<_T, _Alloc>::Matrix(const AbstractMatrix<_U, _Alloc2> &copy) {
	printf("Copied!\n");
	this->dimensions = copy.getDimension();
	this->shape = new int[this->dimensions];
	for (int i = 0; i < this->dimensions; i++) {
		this->shape[i] = copy.getShape(i);
	}
	this->size = copy.getSize();
	this->data = this->allocator.allocate(this->size);
	std::vector<int> index(this->dimensions);
	for (int i = 0; i < this->size; i++) {
		int hold = i;
		//TODO optimize this.
		for(int j = this->dimensions - 1; j >= 0; j--) {
			index[j] = hold % this->getShape(j);
			hold /= this->getShape(j);
		}
		this->data[i] = (_T) copy.getEntry(index);
	}
}

template<typename _T, typename _Alloc>
compchem::Matrix<_T, _Alloc>::Matrix(_T *data,
		std::initializer_list<int> shape) {
	this->dimensions = shape.size();
	this->shape = new int[this->dimensions];
	const int *list = shape.begin();
	for (int i = 0; i < this->dimensions; i++) {
		this->shape[i] = list[i];
	}
	int prod = 1;
	for (int i = 0; i < dimensions; i++) {
		prod *= this->shape[i];
	}
	this->size = prod;
	this->data = this->allocator.allocate(this->size);
	for (int i = 0; i < prod; i++) {
		this->data[i] = data[i];
	}
}

template<typename _T, typename _Alloc>
compchem::Matrix<_T, _Alloc>::~Matrix() {
	this->allocator.deallocate(this->data, this->size);
	delete[] shape;
}

template<typename _T, typename _Alloc>
const _T &compchem::Matrix<_T, _Alloc>::getEntry(
		int index, ...) const {
	int offset = index;
	va_list list;
	va_start(list, index);

	for (int i = 1; i < this->dimensions; i++) {
		int ind = va_arg(list, int);
		if (ind >= this->shape[i] || ind < 0) {
			throw(std::out_of_range(
					"Error: Index out of range. Expected between 0 and "
					+ std::to_string(this->shape[i]) + ", got "
					+ std::to_string(ind) + "\n"));
		}
		offset *= this->getShape(i);
		offset += ind;
	}
	va_end(list);
	return (this->data[offset]);
}

template<typename _T, typename _Alloc>
const _T &compchem::Matrix<_T, _Alloc>::getEntry(
		std::vector<int> index) const {
	if (this->dimensions != index.size()) {
		throw(std::length_error(
				"Error: Attempted to index with incorrect dimensions. Expected: "
						+ std::to_string(this->dimensions) + "; got: "
						+ std::to_string(index.size()) + "\n"));
	}
	int offset = 0;
	int prod = 1;

	for (int i = this->dimensions - 1; i >= 0; i--) {
		if (index[i] >= this->shape[i] || this->shape[i] < 0) {
			throw(std::out_of_range(
					"Error: Index out of range. Expected between 0 and "
							+ std::to_string(this->shape[i]) + ", got "
							+ std::to_string(index[i]) + "\n"));
		}
		offset += index[i] * prod;
		prod *= this->shape[i];
	}
	return (this->data[offset]);
}

template<typename _T, typename _Alloc>
void compchem::Matrix<_T, _Alloc>::setEntry(const _T &ent,
		int index, ...) {
	int offset = index;
	va_list list;
	va_start(list, index);

	for (int i = 1; i < this->dimensions; i++) {
		int ind = va_arg(list, int);
		if (ind >= this->shape[i] || ind < 0) {
			throw(std::out_of_range(
					"Error: Index out of range. Expected between 0 and "
							+ std::to_string(this->shape[i]) + ", got "
							+ std::to_string(index) + "\n"));
		}
		offset *= this->getShape(i);
		offset += ind;
	}
	this->data[offset] = ent;
	va_end(list);
}

template<typename _T, typename _Alloc>
void compchem::Matrix<_T, _Alloc>::setEntry(const _T &ent,
		std::vector<int> index) {
	if (this->dimensions != index.size()) {
		throw(std::length_error(
				"Error: Attempted to index with incorrect dimensions. Expected: "
						+ std::to_string(this->dimensions) + "; got: "
						+ std::to_string(index.size()) + "\n"));
	}
	int offset = 0;
	int prod = 1;

	for (int i = this->dimensions - 1; i >= 0; i--) {
		if (index[i] >= this->shape[i] || this->shape[i] < 0) {
			throw(std::out_of_range(
					"Error: Index out of range. Expected between 0 and "
							+ std::to_string(this->shape[i]) + ", got "
							+ std::to_string(index[i]) + "\n"));
		}
		offset += index[i] * prod;
		prod *= this->shape[i];
	}
	this->data[offset] = ent;
}

template<typename _T, typename _Alloc>
void compchem::Matrix<_T, _Alloc>::setEntry(const _T &&ent,
		int index, ...) {
	int offset = index;
	va_list list;
	va_start(list, index);

	for (int i = 1; i < this->dimensions; i++) {
		int ind = va_arg(list, int);
		if (ind >= this->shape[i] || ind < 0) {
			throw(std::out_of_range(
					"Error: Index out of range. Expected between 0 and "
							+ std::to_string(this->shape[i]) + ", got "
							+ std::to_string(index) + "\n"));
		}
		offset *= this->getShape(i);
		offset += ind;
	}
	this->data[offset] = ent;
	va_end(list);
}

template<typename _T, typename _Alloc>
void compchem::Matrix<_T, _Alloc>::setEntry(const _T &&ent,
		std::vector<int> index) {
	if (this->dimensions != index.size()) {
		throw(std::length_error(
				"Error: Attempted to index with incorrect dimensions. Expected: "
						+ std::to_string(this->dimensions) + "; got: "
						+ std::to_string(index.size()) + "\n"));
	}
	int offset = 0;
	int prod = 1;

	for (int i = this->dimensions - 1; i >= 0; i--) {
		if (index[i] >= this->shape[i] || this->shape[i] < 0) {
			throw(std::out_of_range(
					"Error: Index out of range. Expected between 0 and "
							+ std::to_string(this->shape[i]) + ", got "
							+ std::to_string(index[i]) + "\n"));
		}
		offset += index[i] * prod;
		prod *= this->shape[i];
	}
	this->data[offset] = ent;
}

template<typename _T>
static int compare(const void *a, const void *b) {
	if(*(_T *) a > *(_T *) b) {
		return (1);
	} else if(*(_T *) a < *(_T *) b) {
		return (-1);
	} else {
		return (0);
	}
}

template<typename _T>
static int rcompare(const void *a, const void *b) {
	if(*(_T *) a > *(_T *) b) {
		return (-1);
	} else if(*(_T *) a < *(_T *) b) {
		return (1);
	} else {
		return (0);
	}
}

template<typename _T, typename _Alloc>
void compchem::Matrix<_T, _Alloc>::sort() {
	qsort(this->data, this->size, sizeof(_T), compare<_T>);
}

template<typename _T, typename _Alloc>
void compchem::Matrix<_T, _Alloc>::rsort() {
	qsort(this->data, this->size, sizeof(_T), rcompare<_T>);
}

#endif
