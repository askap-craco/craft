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
#include <stdarg.h>

#define MAX_DIMS 5

typedef float fdmt_dtype;

typedef struct _array2df
{
  int nx;
  int ny;
  fdmt_dtype* d;
} array2d_t;

typedef struct _array3df
{
  int nx;
  int ny;
  int nz;
  fdmt_dtype* d;
} array3d_t;

typedef struct _array4df
{
  int nw;
  int nx;
  int ny;
  int nz;
  fdmt_dtype* d;
} array4d_t;

typedef struct _arraynd
{
	int shape[MAX_DIMS];
	int ndim;
	fdmt_dtype* d;
} arraynd_t;

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

int array4d_idx(const array4d_t* a, int w, int x, int y, int z)
{
  assert(w >=0 && w < a->nw);
  assert(x >=0 && x < a->nx);
  assert(y >=0 && y < a->ny);
  assert(z >=0 && z < a->nz);
  int idx = z + a->nz*(y + a->ny*(x + w*a->nx));
  return idx;
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

#endif /* ARRAY_H_ */
