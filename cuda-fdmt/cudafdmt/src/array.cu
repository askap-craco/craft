/*
 * array.c
 *
 *  Created on: 4 Oct 2016
 *      Author: ban115
 */

#include "array.h"
#include <stdarg.h>


int arraynd_idx(const arraynd_t *a, ...)
{
	va_list argp;
	va_start(argp, a);
	assert(a->ndim >= 1);
	assert(a->ndim <= MAX_DIMS);
	int idx = 0;
	for (int d = 1; d <= a->ndim; d++) {
		int i = va_arg(argp, int);
		int dimsize = a->shape[d];
		assert(i >= 0);
		assert(i < dimsize);
		if (d < a->ndim) {
			int next_dim_size = a->shape[d+1];
			idx = i + next_dim_size*idx;
		} else {
			idx += i;
		}
	}

	va_end(argp);
	return idx;
}

__host__ __device__ int array4d_idx(const array4d_t* a, int w, int x, int y, int z)
{
  //assert(w >=0 && w < a->nw);
//  assert(x >=0 && x < a->nx);
//  assert(y >=0 && y < a->ny);
//  assert(z >=0 && z < a->nz);
  int idx = z + a->nz*(y + a->ny*(x + w*a->nx));
  return idx;
}

__host__ __device__ int array4d_idx(int nw, int nx, int ny, int nz, int w, int x, int y, int z)
{
	  int idx = z + nz*(y + ny*(x + w*nx));
	  return idx;
}

__host__ __device__ size_t array4d_size(const array4d_t* a)
{
	return a->nw * a->nx * a->ny * a->nz;
}

__host__ __device__ size_t array2d_size(const array2d_t* a)
{
	return  a->nx * a->ny;
}

size_t array4d_malloc_hostonly(array4d_t* a)
{
	size_t size = array4d_size(a);
	a->d = (fdmt_dtype*) malloc(size*sizeof(fdmt_dtype));
	assert(a->d != NULL);
	return size;
}

size_t array4d_malloc(array4d_t* a, bool host, bool device)
{
	size_t size = 0;
	if (host) {
		size = array4d_malloc_hostonly(a);
	} else {
		a->d = NULL;
	}

	if (device) {
		size  = array4d_size(a);
		size_t free, total;
		gpuErrchk(cudaMemGetInfo(&free, &total));

		printf("Allocating [%d, %d, %d, %d] %d MIB total %d/%d on GPU\n",a->nw, a->nx, a->ny, a->nz,
				size*sizeof(fdmt_dtype)/1024/1024, (total-free)/1024/1024, total/1024/1024);

		gpuErrchk( cudaMalloc((void**) &a->d_device,
				size*sizeof(fdmt_dtype) ));
		gpuErrchk(cudaMemGetInfo(&free, &total));

		printf("Allocated [%d, %d, %d, %d] %d MIB total %d/%d on GPU\n",a->nw, a->nx, a->ny, a->nz,
		size*sizeof(fdmt_dtype)/1024/1024, (total-free)/1024/1024, total/1024/1024);
	} else {
		a->d_device = NULL;
	}
    return size;
}



int array2d_malloc_hostonly(array2d_t* a)
{
	int size = array2d_size(a);
	a->d = (fdmt_dtype*) malloc(size*sizeof(fdmt_dtype));
    assert(a->d != NULL);
    gpuErrchk( cudaMalloc((void**) &a->d_device, size*sizeof(fdmt_dtype) ));
    return size;
}

int array2d_malloc(array2d_t* a)
{
	int size = array2d_size(a);
	a->d = (fdmt_dtype*) malloc(size*sizeof(fdmt_dtype));
    assert(a->d != NULL);
    gpuErrchk( cudaMalloc((void**) &a->d_device, size*sizeof(fdmt_dtype) ));
    return size;
}

int array4d_copy_to_host(array4d_t* a)
{
	size_t size = array4d_size(a);
	assert(a->d != NULL);
	assert(a->d_device != NULL);
	gpuErrchk(cudaMemcpy(a->d, a->d_device, size*sizeof(fdmt_dtype), cudaMemcpyDeviceToHost));
	return size;
}

int array4d_cuda_memset(array4d_t*a, char c) {
	size_t size = array4d_size(a);
	gpuErrchk(cudaMemset(a->d_device, c, size*sizeof(fdmt_dtype)));
	return size;
}

int array2d_copy_to_device(array2d_t* a)
{
	size_t size = array2d_size(a);
	gpuErrchk(cudaMemcpy(a->d_device, a->d, size*sizeof(fdmt_dtype), cudaMemcpyHostToDevice));
	return size;
}

int array4d_copy_to_device(array4d_t* a)
{
	size_t size = array4d_size(a);
	assert(a->d_device != NULL);
	assert(a->d != NULL);
	gpuErrchk(cudaMemcpy(a->d_device, a->d, size*sizeof(fdmt_dtype), cudaMemcpyHostToDevice));
	return size;
}

int array3d_idx(const array3d_t* a, int x, int y, int z)
{
  assert(x >=0 && x < a->nx);
  assert(y >=0 && y < a->ny);
  assert(z >=0 && z < a->nz);
  int idx = z + a->nz*(y + a->ny*x);
  return idx;
}

int array2d_idx(const array2d_t* a, int x, int y)
{
  assert(x >=0 && x < a->nx);
  assert(y >=0 && y < a->ny);
  if (!(y >=0 && y < a->ny)) {
  }

  int idx = y + a->ny*x;

  return idx;
}

int array2d_dump(const array2d_t* a, const char* foutname)
{
  FILE* fout = fopen(foutname, "w");
  fwrite(&a->nx, sizeof(int), 1, fout);
  fwrite(&a->ny, sizeof(int), 1, fout);
  fwrite(a->d, sizeof(fdmt_dtype), a->nx*a->ny, fout);
  fclose(fout);

  return 0;
}

int array3d_dump(const array3d_t* a, const char* foutname)
{
  FILE* fout = fopen(foutname, "w");
  fwrite(&a->nx, sizeof(int), 1, fout);
  fwrite(&a->ny, sizeof(int), 1, fout);
  fwrite(&a->nz, sizeof(int), 1, fout);
  fwrite(a->d, sizeof(fdmt_dtype), a->nx*a->ny*a->nz, fout);
  fclose(fout);
  return 0;
}

void array4d_print_shape(const array4d_t* a)
{
	printf("nw=%d nx=%d ny=%d nz=%d\b", a->nw, a->nx, a->ny, a->nz);
}

void array4d_set(array4d_t* a, fdmt_dtype v)
{
	assert(a->d != NULL);
	for(int i = 0; i < array4d_size(a); ++i) {
		a->d[i] = v;
	}
	array4d_copy_to_device(a);
}

size_t array4d_zero(array4d_t* a) {
	size_t size = array4d_size(a);
	if (a->d_device) {
		gpuErrchk(cudaMemset(a->d_device, size*sizeof(fdmt_dtype), 0));
	}
	if (a->d) {
		bzero(a->d, size*sizeof(fdmt_dtype));
	}

	return size;
}

int array4d_dump(const array4d_t* a, const char* foutname)
{
  FILE* fout = fopen(foutname, "w");
  fwrite(&a->nw, sizeof(int), 1, fout);
  fwrite(&a->nx, sizeof(int), 1, fout);
  fwrite(&a->ny, sizeof(int), 1, fout);
  fwrite(&a->nz, sizeof(int), 1, fout);
  fwrite(a->d, sizeof(fdmt_dtype), a->nw*a->nx*a->ny*a->nz, fout);
  fclose(fout);
  return 0;
}
