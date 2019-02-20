/*
 * disk_allocator.hpp
 *
 *  Created on: Feb 18, 2019
 *      Author: connor
 */

#ifndef PROJECT_11_DISK_ALLOCATOR_HPP_
#define PROJECT_11_DISK_ALLOCATOR_HPP_

#include <stdlib.h>
#include <map>
#include <tuple>

namespace compchem {

namespace strategies {

template<typename _T>
class DiskAllocator {
public:
	typedef _T value_type;
	typedef _T *pointer;

	DiskAllocator();

	virtual ~DiskAllocator();

	pointer allocate(size_t n_elems);

	void deallocate(pointer ptr, size_t n_elems);
private:
#ifdef __WIN32
	std::map<const pointer, size_t> regions;
#else
	std::map<const pointer, std::pair<int, size_t> > files;
#endif
};

}

}

#include "disk_allocator.cpp"

#endif /* PROJECT_11_DISK_ALLOCATOR_HPP_ */
