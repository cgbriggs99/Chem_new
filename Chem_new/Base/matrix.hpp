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

namespace compchem {

template<typename _T, typename _Alloc = std::allocator<_T> >
class Matrix {
public:
  //Constructors and destructor.
  Matrix(std::initializer_list<int> shape);
  Matrix(const Matrix<_T> &copy);

  template<typename _U, typename _Alloc2 = std::allocator<_U> >
  Matrix(const Matrix<_U, _Alloc2> &copy_and_cast);
  Matrix(_T *data, std::initializer_list<int> shape);

  virtual ~Matrix();

protected:
  //Data
  _T *data;
  std::vector<int> shape;
  int dimensions;
  int size;
  _Alloc allocator;
public:
  //Getters and setters.
  virtual _T &getEntry(std::initializer_list<int> index);
  virtual const _T &getEntry(std::initializer_list<int> index) const;
  virtual int getSize() {
    return (this->size);
  }
  
  virtual const int *getShape() {
    return (this->shape);
  }

  virtual int getDimension() {
    return (this->dimensions);
  }
};
}

#include "matrix.cpp"

#endif
