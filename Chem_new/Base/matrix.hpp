/*
 * matrix.hpp
 *
 * Created on: Dec 22, 2018
 *     Author: Connor
 */

#ifndef MATRIX_HPP_
#define MATRIX_HPP_

#include <memory>
#include <initializer_list>
#include <vector>
#include <stdarg.h>

namespace compchem {

template<typename T, typename _Alloc = std::allocator<T> >
class AbstractMatrix {
public:
	AbstractMatrix() {
		;
	}
	virtual ~AbstractMatrix() = default;

	virtual const T &getEntry(int index, ...) const = 0;
	virtual const T &getEntry(std::vector<int> index) const = 0;
	virtual void setEntry(const T &ent, int index, ...) = 0;
	virtual void setEntry(const T &ent, std::vector<int> index) = 0;
	virtual void setEntry(const T &&ent, int index, ...) = 0;
	virtual void setEntry(const T &&ent, std::vector<int> index) = 0;
	virtual int getSize() const = 0;
	virtual int getShape(int dim) const = 0;
	virtual int getDimension() const = 0;
};

}

#endif
