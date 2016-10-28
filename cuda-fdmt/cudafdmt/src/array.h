/*
 * array.h
 *
 *  Created on: 22 Sep 2016
 *      Author: ban115
 */

#ifndef ARRAY_H_
#define ARRAY_H_

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include "cuda_utils.h"

#define MAX_DIMS 5

typedef float fdmt_dtype;

typedef struct _array2df
{
  int nx;
  int ny;
  fdmt_dtype* d;
  fdmt_dtype* d_device;
} array2d_t;

typedef struct _array3df
{
  int nx;
  int ny;
  int nz;
  fdmt_dtype* d;
  fdmt_dtype* d_device;
} array3d_t;

typedef struct _array4df
{
  int nw;
  int nx;
  int ny;
  int nz;
  fdmt_dtype* d;
  fdmt_dtype* d_device;
} array4d_t;

typedef struct _arraynd
{
	int shape[MAX_DIMS];
	int ndim;
	fdmt_dtype* d;
	fdmt_dtype* d_device;
} arraynd_t;

typedef struct _coord3 {
	int x;
	int y;
	int z;
} coord3_t;

typedef struct _coord4 {
	int w;
	int x;
	int y;
	int z;
} coord4_t;

int arraynd_idx(const arraynd_t *a, ...);

size_t array4d_malloc(array4d_t* a);
size_t array4d_malloc_hostonly(array4d_t* a);
int array4d_copy_to_host(array4d_t* a);
int array4d_cuda_memset(array4d_t*a, char c) ;
int array4d_copy_to_device(array4d_t* a);
void array4d_print_shape(array4d_t* a);


int array2d_malloc(array2d_t* a);
int array2d_copy_to_device(array2d_t* a);

__host__ __device__ size_t array4d_size(const array4d_t* a);
__host__ __device__ int array4d_idx(const array4d_t* a, int w, int x, int y, int z);
__host__ __device__ int array3d_idx(const array3d_t* a, int x, int y, int z);
__host__ __device__ int array2d_idx(const array2d_t* a, int x, int y);

int array2d_dump(const array2d_t* a, const char* foutname);
int array3d_dump(const array3d_t* a, const char* foutname);
int array4d_dump(const array4d_t* a, const char* foutname);

#endif /* ARRAY_H_ */
