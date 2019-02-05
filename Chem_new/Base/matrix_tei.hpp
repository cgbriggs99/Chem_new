/*
 * matrix_tei.hpp
 *
 *  Created on: Feb 4, 2019
 *      Author: Connor
 */

#ifndef BASE_MATRIX_TEI_HPP_
#define BASE_MATRIX_TEI_HPP_

#include <memory>
#include "matrix.hpp"

namespace compchem {
namespace strategies {

template<typename _T, typename _Alloc = std::allocator<_T> >
class TEIMatrix : public AbstractMatrix<_T, _Alloc> {
private:
	_T *data;
	int n;
	int size;

public:
	TEIMatrix(int n);
	virtual ~TEIMatrix();

	const _T &getEntry(int index, ...) const;
	const _T &getEntry(std::vector<int> index) const;
	void setEntry(const _T &ent, int index, ...);
	void setEntry(const _T &ent, std::vector<int> index);
	void setEntry(const _T &&ent, int index, ...);
	void setEntry(const _T &&ent, std::vector<int> index);
	int getSize() const;
	int getShape(int dim) const;
	int getDimension() const;
};

}
}

#include "matrix_tei.cpp"

#endif /* BASE_MATRIX_TEI_HPP_ */
