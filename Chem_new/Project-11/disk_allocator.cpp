/*
 * disk_allocator.cpp
 *
 *  Created on: Feb 18, 2019
 *      Author: connor
 */

#ifndef __DISK_ALLOCATOR_CPP__
#define __DISK_ALLOCATOR_CPP__

#include "disk_allocator.hpp"

#include <stdlib.h>

#ifdef __WIN32
#include <windows.h>

template<typename _T>
compchem::strategies::DiskAllocator<_T>::DiskAllocator() {

}

template<typename _T>
compchem::strategies::DiskAllocator<_T>::~DiskAllocator() {
	for(std::pair<pointer, size_t> &pair : this->regions) {
		VirtualFree(pair.first, pair.second, MEM_RELEASE);
	}
}

template<typename _T>
compchem::strategies::DiskAllocator<_T>::pointer compchem::strategies::DiskAllocator<_T>::allocate(size_t n_elems) {
	size_t size = n_elems * sizeof(_T);
	pointer ptr = VirtualAlloc(NULL, size, MEM_RESERVE | MEM_COMMIT, PAGE_READWRITE);
	regions[ptr] = size;
}

template<typename _T>
void compchem::strategies::DiskAllocator<_T>::deallocate(pointer ptr, size_t n_elems) {
	if(regions.find(ptr) == regions.end()) {
		return;
	}

	if(n_elems * sizeof(_T) != regions[ptr].second) {
		return;
	}
	VirtualFree(regions[ptr].first, regions[ptr].second);
	regions.erase(ptr);
}

#else
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>

template<typename _T>
compchem::strategies::DiskAllocator<_T>::DiskAllocator() {

}

template<typename _T>
compchem::strategies::DiskAllocator<_T>::~DiskAllocator() {
	for(std::pair<const pointer, std::pair<int, size_t> > &pair : this->files) {
		munmap(pair.first, pair.second.second);
		close(pair.second.first);
	}
}

#define PAGE_SIZE 4196

template<typename _T>
typename compchem::strategies::DiskAllocator<_T>::pointer compchem::strategies::DiskAllocator<
        _T>::allocate(size_t n_elems) {
	size_t size = n_elems * sizeof(_T);

	//Open the file in a temporary directory.
	int fd = open("/tmp", O_TMPFILE | O_RDWR);

	//Fill the file initially.
	char buff[PAGE_SIZE] = {0};
	for(int i = 0; i < size / PAGE_SIZE; i++) {
		write(fd, buff, PAGE_SIZE);
	}
	if(size % PAGE_SIZE != 0) {
		write(fd, buff, size % PAGE_SIZE);
	}

	//Map the file into virtual memory.
	pointer ptr = (pointer) mmap(NULL, size, PROT_READ | PROT_WRITE,
	        MAP_SHARED | MAP_POPULATE, fd, 0);
	//Add file details to the map.
	files[ptr] = std::pair<int, size_t>(fd, size);
	return (ptr);
}

template<typename _T>
void compchem::strategies::DiskAllocator<_T>::deallocate(pointer ptr,
        size_t n_elems) {
	if(files.find(ptr) == files.end()) {
		return;
	}

	if(n_elems * sizeof(_T) != files[ptr].second) {
		return;
	}

	munmap(ptr, files[ptr].second);
	close(files[ptr].first);
	files.erase(ptr);
}

#undef PAGE_SIZE

#endif

#endif
