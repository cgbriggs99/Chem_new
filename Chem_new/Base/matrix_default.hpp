/*
 * matrix_default.hpp
 *
 *  Created on: Feb 4, 2019
 *      Author: connor
 */

#ifndef BASE_MATRIX_DEFAULT_HPP_
#define BASE_MATRIX_DEFAULT_HPP_

#include <memory>
#include <initializer_list>
#include <vector>
#include "matrix.hpp"

namespace compchem {

template<typename _T, typename _Alloc = std::allocator<_T> >
class Matrix: public AbstractMatrix<_T, _Alloc> {
public:
	//Constructors and destructor.
	Matrix(std::initializer_list<int> shape);
	explicit Matrix(const AbstractMatrix<_T> &copy);

	template<typename _U, typename _Alloc2 = std::allocator<_U> >
	explicit Matrix(const AbstractMatrix<_U, _Alloc2> &copy_and_cast);
	Matrix(_T *data, std::initializer_list<int> shape);

	virtual ~Matrix();

protected:
	//Data
	_T *data;
	int *shape;
	int dimensions;
	int size;
	_Alloc allocator;
public:
	//Getters and setters.
	const _T &getEntry(int index, ...) const override;
	const _T &getEntry(std::vector<int> index) const override;
	void setEntry(_T &ent, int index, ...) override;
	void setEntry(_T &ent, std::vector<int> index) override;
	void setEntry(_T &&ent, int index, ...) override;
	void setEntry(_T &&ent, std::vector<int> index) override;
	int getSize() const override {
		return (this->size);
	}

	int getShape(int ind) const override {
		return (this->shape[ind]);
	}

	int getDimension() const override {
		return (this->dimensions);
	}

	//Sorts low to high.
	void sort();
	//Sorts high to low.
	void rsort();
};
}

#include "matrix.cpp"

#endif /* BASE_MATRIX_DEFAULT_HPP_ */
