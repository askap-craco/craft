/*
 * fftbench.cu
 *
 *  Created on: 17 Sep 2018
 *      Author: ban115
 *  Modified by Xinping Deng at August 2019
 *
 */


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cuda.h>
#include <cufft.h>
#include <cufftXt.h>
#include <cuda_fp16.h>
#include "CudaTimer.h"
#include "CpuTimer.h"
#include "cuda_utils.h"
#include "cufft_utils.h"

__device__ cufftComplex load_callback(void *dataIn,
				      size_t offset,
				      void *callerInfo,
				      void *sharedPtr) {
  cufftComplex value = {0.0f, 0.0f};
  return value;
}
__device__ cufftCallbackLoadC d_loadCallbackPtr = load_callback;

__device__ void store_callback(void *dataOut,
			       size_t offset,
			       cufftReal element,
			       void *callerInfo,
			       void *sharedPtr) {
}
__device__ cufftCallbackStoreR d_storeCallbackPtr = store_callback;

template <class intype>
int timefft(int n, int batch, cudaDataType itype, cudaDataType etype, cudaDataType otype, int option, bool inplace)
{
  CudaTimer t;
  intype *data, *out_data, *data_host;
  cufftHandle plan;
  size_t data_size=sizeof(intype)*n*(n/2 + 1)*batch;

  gpuErrchk(cudaMalloc((void**) &data, data_size));
  if (inplace) {
    out_data = data;
  } else {
    gpuErrchk(cudaMalloc((void**) &out_data, data_size));
  }
  size_t worksize;
  cufftSafeCall(cufftCreate(&plan));
  
  /* To get numbers in the array */
  srand(time(NULL));
  gpuErrchk(cudaMallocHost((void **)&data_host, data_size));
  for(int i = 0; i < n*(n/2+1); i++)
    {
      data_host[i].x = (float)rand();
      data_host[i].y = (float)rand();
    }
  gpuErrchk(cudaMemcpy(data, data_host, data_size, cudaMemcpyHostToDevice));
  
  if(option != 1)
    {
      long long int nsize[] = {n,n };
      
      cufftSafeCall(cufftXtMakePlanMany(plan, 2, nsize,
					NULL, 1, 0, itype,
					NULL, 1, 0, otype,
					batch, &worksize, etype
					));
    }
  else
    {
      long long int nsize = n;
      cufftSafeCall(cufftXtMakePlanMany(plan, 1, &nsize,
					NULL, 1, 0, itype,
					NULL, 1, 0, otype,
					batch * n, &worksize, etype
					));
    }

  if(((option == 2) || (option == 3) || (option == 4)) && (sizeof(intype) != 8))
    {
      printf("%dx%d FFT batch=%d data=%d MB in-place=%d type=%d-> %d. cuFFT callback only support single precision\n",
	     n,n,batch,data_size/1024/1024, inplace, itype,otype);
      
      cufftSafeCall(cufftDestroy(plan));
      gpuErrchk(cudaFree(data));
      gpuErrchk(cudaFreeHost(data_host));
      if (! inplace) {
	gpuErrchk(cudaFree(out_data));
      }
      return EXIT_SUCCESS;
    }
  
  if((option == 2) || (option == 4))
    {
      cufftCallbackLoadR h_loadCallbackPtr;
      gpuErrchk(cudaMemcpyFromSymbol(&h_loadCallbackPtr,
  				 d_loadCallbackPtr, 
  				 sizeof(h_loadCallbackPtr)));
      // Now associate the callbacks with the plan.
      cufftResult status = cufftXtSetCallback(plan, 
					      (void **)&h_loadCallbackPtr, 
					      CUFFT_CB_LD_COMPLEX,
					      0);
      if (status == CUFFT_LICENSE_ERROR) {
	fprintf(stdout, "This sample requires a valid license file.\n");
	fprintf(stdout, "The file was either not found, out of date, or otherwise invalid.\n");
	exit(EXIT_FAILURE);
      } else {
	cufftSafeCall(status);
      }
    }
  if((option == 3) || (option == 4))
    {
      cufftCallbackStoreR h_storeCallbackPtr;
      gpuErrchk(cudaMemcpyFromSymbol(&h_storeCallbackPtr,
				 d_storeCallbackPtr,
				 sizeof(h_storeCallbackPtr)));
      cufftSafeCall(cufftXtSetCallback(plan,
				       (void **)&h_storeCallbackPtr,
				       CUFFT_CB_ST_REAL,
				       NULL));
    }

  // warm up
  cufftSafeCall(cufftXtExec(plan, data, out_data, CUFFT_INVERSE));

  int niter = 100;
  for (int i = 0; i < niter; ++i) {
    t.start();
    cufftSafeCall(cufftXtExec(plan, data, out_data, CUFFT_INVERSE));
    t.stop();
  }
  
  float tavg_us = t.get_average_time() / float(batch) * 1e3f;

  printf("%dx%d FFT batch=%d data=%d MB in-place=%d type=%d-> %d. Worksize=%d MB: %f microseconds/FFT= %f k FFTs/sec total=%0.2fs\n",
	 n,n,batch,data_size/1024/1024, inplace, itype,otype, worksize/1024/1024, tavg_us, 1./tavg_us*1e6f/1e3f);
  cufftSafeCall(cufftDestroy(plan));
  gpuErrchk(cudaFree(data));
  gpuErrchk(cudaFreeHost(data_host));
  if (! inplace) {
    gpuErrchk(cudaFree(out_data));
  }

  return EXIT_SUCCESS;
}

int main(int argc, char* argv[])
{

  if (argc != 6) {
    printf("%s Usage: gpuid N batchmin batchmax option\n", argv[0]);
    return EXIT_FAILURE;
  }

  int n = atoi(argv[2]);
  int batchmin = atoi(argv[3]);
  int batchmax = atoi(argv[4]);
  int option   = atoi(argv[5]);

  int cuda_device = atoi(argv[1]);
  cudaDeviceProp p;
  gpuErrchk(cudaGetDeviceProperties(&p, cuda_device));

  if(option == 0)
    printf("Benchmark with original 2D FFT \n");
  else if(option == 1)
    printf("Benchmark with 1D FFT \n");
  else if(option == 2)
    printf("Benchmark with 2D FFT using load callback \n");
  else if(option == 3)
    printf("Benchmark with 2D FFT using store callback \n");
  else if(option == 4)
    printf("Benchmark with 2D FFT using load and store callback \n");
  else
    {
      printf("The option can only be 0 to 4\n");
      return EXIT_FAILURE;
    }
    
  printf("Device[%d]=%s v%d.%d Mem=%d GB shmem/block=%d constmem=%d Warp=%d Clock=%d MHz %d multiprocessors\n",
	 cuda_device, p.name, p.major, p.minor, p.totalGlobalMem/1024/1024/1024, p.sharedMemPerBlock,
	 p.totalConstMem, p.warpSize, p.clockRate/1000, p.multiProcessorCount
	 );

  gpuErrchk( cudaSetDevice(cuda_device));

  /* 32 bits complex */
  cudaDataType itype = CUDA_C_32F;
  cudaDataType etype = CUDA_C_32F;
  cudaDataType otype = CUDA_R_32F;

  for (int batch2 = batchmin; batch2 < batchmax; batch2++) {
    int batch = 1 << batch2;
    timefft<cufftComplex>(n,batch,itype, etype, otype, option, false);
  }
  for (int batch2 = batchmin; batch2 < batchmax; batch2++) {
    int batch = 1 << batch2;
    timefft<cufftComplex>(n,batch,itype, etype, otype, option, true);
  }

  /* Half complex */
  itype = CUDA_C_16F;
  etype = CUDA_C_16F;
  otype = CUDA_R_16F;

  for (int batch2 = batchmin; batch2 < batchmax; batch2++) {
    int batch = 1 << batch2;
    timefft<half2>(n,batch,itype, etype, otype, option, false);
  }  
  for (int batch2 = batchmin; batch2 < batchmax; batch2++) {
    int batch = 1 << batch2;
    timefft<half2>(n,batch,itype, etype, otype, option, true);
  }

  printf("Benchmark finished\n");
}
