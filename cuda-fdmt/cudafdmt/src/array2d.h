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

	Array1d(size_t nx = 0);
	virtual ~Array1d();
	virtual __host__ size_t copy_to_device();
	virtual __host__ size_t copy_to_host();
	__host__ __device__ size_t size() {
		return m_size;
	}
	virtual T &operator[](int x) {
		return d[x];
	}
};

template<class T>
class Array2d : public Array1d<T>
{
	int m_nx;
	int m_ny;

public:
	Array2d();
	Array2d(int nx, int ny);

	__host__ int index(int x, int y) const {
		assert(x >= 0);
		assert(y >= 0);
		assert(x < m_nx);
		assert(y < m_ny);
		return m_ny*x + y;
	}

	__host__ __device__ T &operator()(int x, int y) {
		return this->d[index(x, y)]; // <--- Compiler is happy
	}

	__host__ void set_host(int x, int y, T& value) {
		this->d[index(x, y)] = value; // <--- Compiler complains '"d" is undefined". Works if I do this->d[index(x, y)] = value;
	}

	size_t dump(const char* foutname) const; // Implemented down the bottom of my file

};

template <class T>
Array1d<T>::Array1d(size_t sz) {
	if (sz == 0) {
		d = NULL;
		d_device = NULL;
		m_size = 0;
	} else {
		printf("Allocating Array1d size %d\n", sz);
		d = new T[sz];
		assert(d != NULL);
		printf("CUDA Allocating Array1d size %d\n", sz);

		gpuErrchk( cudaMalloc((void**) &d_device, sz*sizeof(T) ));
		m_size = sz;
		printf("Alocated Array1d size %d\n", sz);

	}

}

template <class T>
Array1d<T>::~Array1d() {
	if (d != NULL) {
		printf("Freeing 0x%x\n", d);
		delete[] d;
		printf("Freeed0x%x\n", d);
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
	return sizeof(T)*this->m_size;
}

template <class T>
__host__ size_t Array1d<T>::copy_to_device() {
	size_t sz = m_size;
	printf("Copying to ddevive %d host address 0x%x device address: 0x%x\n", sz, d, d_device);
	gpuErrchk(cudaMemcpy(d_device, d, sz*sizeof(T), cudaMemcpyHostToDevice));
	return sz*sizeof(T);
}

template <class T>
__host__ size_t Array1d<T>::copy_to_host() {
	size_t sz = m_size;
	gpuErrchk(cudaMemcpy(d, d_device, sz*sizeof(T), cudaMemcpyDeviceToHost));
	return sz;
}

template <class T>
Array2d<T>::Array2d(int nx, int ny ) : Array1d<T>(nx*ny), m_nx(nx), m_ny(ny) {
	printf("Array2d %d %d\n", nx, ny);
}


#endif /* ARRAY2D_H_ */
