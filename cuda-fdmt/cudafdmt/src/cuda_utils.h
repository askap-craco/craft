/*
 * cuda_utils.h
 *
 *  Created on: 27 Sep 2016
 *      Author: ban115
 */

#ifndef CUDA_UTILS_H_
#define CUDA_UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <assert.h>

#ifdef __cplusplus
extern "C"
#endif

__host__  inline void gpuAssert(cudaError_t code, const char *file, int line)
{
   if (code != cudaSuccess)
   {
      fprintf(stderr,"GPUassert: %s %s:%d\n", cudaGetErrorString(code), file, line);
      assert(code == cudaSuccess);
      exit(code);
   }
}


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }



#endif /* CUDA_UTILS_H_ */
