#ifndef __MATRIX_CPP__
#define __MATRIX_CPP__

#include "matrix.hpp"
#include <exception>
#include <stdexcept>
#include <string>
#include <initializer_list>

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
compchem::Matrix<_T, _Alloc>::Matrix(const Matrix<_T> &copy) {
	puts("Copy constructor!");
	this->dimensions = copy.dimensions;
	this->size = copy.size;
	this->shape = new int[this->dimensions];
	for (int i = 0; i < this->dimensions; i++) {
		this->shape[i] = copy.shape[i];
	}
	this->data = this->allocator.allocate(this->size);
	for (int i = 0; i < this->size; i++) {
		this->data[i] = copy.data[i];
	}
}

template<typename _T, typename _Alloc>
template<typename _U, typename _Alloc2>
compchem::Matrix<_T, _Alloc>::Matrix(const Matrix<_U, _Alloc2> &copy) {
	puts("Cast and copy constructor!");
	this->dimensions = copy.dimensions;
	this->shape = new int[this->dimensions];
	for (int i = 0; i < this->dimensions; i++) {
		this->shape[i] = copy.shape[i];
	}
	this->size = copy.size;
	this->data = this->allocator.allocate(this->size);
	for (int i = 0; i < this->size; i++) {
		this->data[i] = (_T) copy.data[i];
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
_T &compchem::Matrix<_T, _Alloc>::getEntry(std::initializer_list<int> index) {
	if (this->dimensions != index.size()) {
		throw(std::length_error(
				"Error: Attempted to index with incorrect dimensions. Expected: "
						+ std::to_string(this->dimensions) + "; got: "
						+ std::to_string(index.size()) + "\n"));
	}
	int offset = 0;
	int prod = 1;
	const int *list = index.begin();
	for (int i = this->dimensions - 1; i >= 0; i--) {
		if (list[i] >= this->shape[i] || this->shape[i] < 0) {
			throw(std::out_of_range(
					"Error: Index out of range. Expected between 0 and "
							+ std::to_string(this->shape[i]) + ", got "
							+ std::to_string(index.begin()[i]) + "\n"));
		}
		offset += list[i] * prod;
		prod *= this->shape[i];
	}
	return (this->data[offset]);
}

template<typename _T, typename _Alloc>
const _T &compchem::Matrix<_T, _Alloc>::getEntry(
		std::initializer_list<int> index) const {
	if (this->dimensions != index.size()) {
		throw(std::length_error(
				"Error: Attempted to index with incorrect dimensions. Expected: "
						+ std::to_string(this->dimensions) + "; got: "
						+ std::to_string(index.size()) + "\n"));
	}
	int offset = 0;
	int prod = 1;
	const int *list = index.begin();

	for (int i = this->dimensions - 1; i >= 0; i--) {
		if (list[i] >= this->shape[i] || this->shape[i] < 0) {
			throw(std::out_of_range(
					"Error: Index out of range. Expected between 0 and "
							+ std::to_string(this->shape[i]) + ", got "
							+ std::to_string(index.begin()[i]) + "\n"));
		}
		offset += list[i] * prod;
		prod *= this->shape[i];
	}
	return (this->data[offset]);
}

#endif
