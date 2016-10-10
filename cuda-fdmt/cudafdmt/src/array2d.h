/*
 * array2d.h
 *
 *  Created on: 10 Oct 2016
 *      Author: ban115
 */

#ifndef ARRAY2D_H_
#define ARRAY2D_H_

#include <cuda.h>
#include <stdio.h>
#include <assert.h>

template<class T> class Array1d
{
public:
	T* d;
	T* d_device;
	size_t m_size;
	Array1d();
	Array1d(size_t nx);
	virtual ~Array1d();
	virtual __host__ size_t copy_to_device() const;
	virtual __host__ size_t copy_to_host();
	__host__ __device__ size_t size() {
		return m_size;
	}
	virtual T &operator[](int x) {
		return d[x];
	}
};

template<class T>
class Array2d : Array1d<T>
{
	int m_nx;
	int m_ny;

public:
	Array2d();
	Array2d(int nx, int ny);
	virtual ~Array2d();
	size_t dump(const char* foutname) const;

	__host__ __device__ T &operator()(int x, int y) {
		return this->d[index(x, y)];
	}
	__host__ __device__ int index(int x, int y) const;

};

template <class T>
Array1d<T>::Array1d() {
	d = NULL;
	d_device = NULL;
	m_size = 0;
}

template <class T>
Array1d<T>::Array1d(size_t sz) {
	d = (T*) malloc(sz*sizeof(T));
	assert(d != NULL);
	gpuErrchk( cudaMalloc((void**) &d_device, sz*sizeof(T) ));
	m_size = sz;
}

template <class T>
Array1d<T>::~Array1d() {
	if (d != NULL) {
		free(d);
	}

	if (d_device != NULL) {
		gpuErrchk( cudaFree(d_device) );
	}
}


template <class T>
__host__ size_t Array2d<T>::dump(const char* foutname) const {
	FILE* fout = fopen(foutname, "w");
	fwrite(&m_nx, sizeof(int), 1, fout);
	fwrite(&m_ny, sizeof(int), 1, fout);
	fwrite(this->d, sizeof(T),  this->m_size , fout);
	fclose(fout);
}

template <class T>
__host__ size_t Array1d<T>::copy_to_device() const {
	size_t sz = m_size;
	gpuErrchk(cudaMemcpy(d_device, d, sz*sizeof(T), cudaMemcpyHostToDevice));
	return sz;
}

template <class T>
__host__ size_t Array1d<T>::copy_to_host() {
	size_t sz = m_size;
	gpuErrchk(cudaMemcpy(d, d_device, sz*sizeof(T), cudaMemcpyDeviceToHost));
	return sz;
}

template <class T>
Array2d<T>::Array2d(int nx, int ny ) : Array1d<T>(nx*ny) {
	m_nx = nx;
	m_ny = ny;
}


template <class T>
__host__ __device__ int Array2d<T>::index(int x, int y) const {
	return m_ny*x + y;
}


#endif /* ARRAY2D_H_ */
