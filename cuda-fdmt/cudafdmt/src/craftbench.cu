/*
 * craftbench.cu
 *
 *  Created on: 17 Sep 2018
 *      Author: ban115
 *  Updated on August 2019 by DEN15C 
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>  
#include <assert.h>
#include <cuda.h>
#include <cufft.h>
#include <cufftXt.h>
#include <cuda_fp16.h>
#include "CudaTimer.h"
#include "CpuTimer.h"
#include "cuda_utils.h"
#include "cufft_utils.h"

void usage()
{
  fprintf(stdout,
	  "craftbench - To check the profermance of CRAFT essentials with different configurations\n"
	  "\n"
	  "Usage: craftbench [options]\n"
	  " -a The ID of GPU to be used \n"
	  " -b FFT size\n"
	  " -c Min number of batch to do, power of 2\n "
	  " -d Max number of batch to do, power of 2\n"
	  " -e Options we want to do with: 0 FFT without host, 1 FFT with memcpy between device and host, more options to be added\n"
	  );
}

/* A kernel to unpack from char to given type */
template <class intype>
__global__ void unpack_kernel(char *dbuf_in, intype *dbuf_out)
{
  uint64_t loc_out = blockIdx.x * gridDim.y * blockDim.x +
    blockIdx.y * blockDim.x +
    threadIdx.x;
  uint64_t loc_in = 2 * loc_out;

  dbuf_out[loc_out].x = dbuf_in[loc_in];
  dbuf_out[loc_out].y = dbuf_in[loc_in + 1];  
}

template <class intype>
int timefft(int n, int batch, bool inplace, int option)
{
  CudaTimer t;
  intype *dbuf_rt, *dbuf_out;
  char *dbuf_in, *hbuf;
  cufftHandle plan;
  size_t dbufrt_size=sizeof(intype)*n*(n/2 + 1)*batch;
  size_t buf_size=2*sizeof(char)*n*(n/2 + 1)*batch;
  cudaDataType itype, etype, otype;
  
  if(sizeof(intype) == sizeof(cufftComplex))
    {
      itype = CUDA_C_32F;
      etype = CUDA_C_32F;
      otype = CUDA_R_32F;
    }
  if(sizeof(intype) == sizeof(half2))
    {
      itype = CUDA_C_16F;
      etype = CUDA_C_16F;
      otype = CUDA_R_16F;
    }
 
  /* Get memory on device */
  gpuErrchk(cudaMalloc((void**) &dbuf_rt, dbufrt_size));
  if (inplace)
    dbuf_out = dbuf_rt;
  else
    gpuErrchk(cudaMalloc((void**) &dbuf_out, dbufrt_size));

  /* Get memory on host */
  if(option != 0)
    {
      gpuErrchk(cudaMalloc((void**) &dbuf_in, buf_size));
      gpuErrchk(cudaMallocHost((void**) &hbuf, buf_size));
    }
  long long int nsize[] = {n,n };

  /* Create the FFT plan for multiple dbuf types */
  size_t worksize;
  cufftSafeCall(cufftCreate(&plan));
  cufftSafeCall(cufftXtMakePlanMany(plan, 2, nsize,
				    NULL, 1, 0, itype,
				    NULL, 1, 0, otype,
				    batch, &worksize, etype
				    ));
  
  /* Warm up and do the real thing */
  int i, niter = 100;
  dim3 gridSize, blockSize;
  blockSize.x = n;
  blockSize.y = 1;
  blockSize.z = 1;
  gridSize.x = n/2 + 1;
  gridSize.y = 1;
  gridSize.z = 1;
  if(option == 1)
    {
      gpuErrchk(cudaMemcpy(dbuf_in, hbuf, buf_size, cudaMemcpyHostToDevice));
      
      unpack_kernel<intype><<<gridSize, blockSize, 0>>>(dbuf_in, dbuf_rt);
    }
  cufftSafeCall(cufftXtExec(plan, dbuf_rt, dbuf_out, CUFFT_INVERSE));
  for (i = 0; i < niter; ++i)
    {
      t.start();
      if(option == 1)
	gpuErrchk(cudaMemcpy(dbuf_rt, hbuf, buf_size, cudaMemcpyHostToDevice));
      cufftSafeCall(cufftXtExec(plan, dbuf_rt, dbuf_out, CUFFT_INVERSE));
      t.stop(); 
    }

  /* Check the timer and display */
  float tavg_us = t.get_average_time() / float(batch) * 1e3f;
  fprintf(stdout, "%dx%d FFT batch=%d dbuf=%d MB in-place=%d type=%d-> %d. Worksize=%d MB: %f microseconds/FFT= %f k FFTs/sec total=%0.2fs\n",
	  n,n,batch,dbufrt_size/1024/1024, inplace, itype,otype, worksize/1024/1024, tavg_us, 1./tavg_us*1e6f/1e3f);

  /* Destroy FFT plan */
  cufftSafeCall(cufftDestroy(plan));

  /* Free device memory */
  gpuErrchk(cudaFree(dbuf_rt));
  if (! inplace)
    gpuErrchk(cudaFree(dbuf_out));

  /* Free host memory */
  if(option != 0)
    {
      gpuErrchk(cudaFreeHost(hbuf));
      gpuErrchk(cudaFree(dbuf_in));
    }
  
  return EXIT_SUCCESS;
}

// ./craftbench -a 0 -b 256 -c 10 -d 12 -e 0
int main(int argc, char* argv[])
{
  int arg;
  int cuda_device, n, batchmax, batchmin, option;
  int narg = 0, narg_expect = 5;
  
  /* read in argument from command line */
  while((arg=getopt(argc,argv,"a:b:c:hd:e:")) != -1)
    {
      switch(arg)
	{
	case 'h':
	  usage();
	  exit(EXIT_FAILURE);
	  
	case 'a':	  	  
	  if(sscanf(optarg, "%d", &cuda_device) != 1)
	    {
	      usage();	      
	      exit(EXIT_FAILURE);
	    }
	  narg++;
	  break;
	  
	case 'b':	  	  
	  if(sscanf(optarg, "%d", &n) != 1)
	    {
	      usage();	      
	      exit(EXIT_FAILURE);
	    }
	  narg++;
	  break;
	  
	case 'c':	  	  
	  if(sscanf(optarg, "%d", &batchmin) != 1)
	    {
	      usage();	      
	      exit(EXIT_FAILURE);
	    }
	  narg++;
	  break;
	  
	case 'd':	  	  
	  if(sscanf(optarg, "%d", &batchmax) != 1)
	    {
	      usage();	      
	      exit(EXIT_FAILURE);
	    }
	  narg++;
	  break;
	  
	case 'e':	  	  
	  if(sscanf(optarg, "%d", &option) != 1)
	    {
	      usage();	      
	      exit(EXIT_FAILURE);
	    }
	  narg++;
	  break;	  
	}
    }

  /* To check if we et all information */
  if (narg != narg_expect)
    {
      usage();	      
      exit(EXIT_FAILURE);
    }

  /* Check available GPU and give information */
  cudaDeviceProp p;
  gpuErrchk(cudaGetDeviceProperties(&p, cuda_device));  
  fprintf(stdout, "FFT Benchmark \n");
  fprintf(stdout, "Device[%d]=%s v%d.%d Mem=%d GB shmem/block=%d constmem=%d Warp=%d Clock=%d MHz %d multiprocessors\n",
	 cuda_device, p.name, p.major, p.minor, p.totalGlobalMem/1024/1024/1024, p.sharedMemPerBlock,
	 p.totalConstMem, p.warpSize, p.clockRate/1000, p.multiProcessorCount
	 );  
  gpuErrchk(cudaSetDevice(cuda_device));

  /* Do the real work */
  int batch2, batch;
  for (batch2 = batchmin; batch2 < batchmax; batch2++)
    {
      batch = 1 << batch2;
      timefft<cufftComplex>(n,batch,false, option);
    }
  
  for (batch2 = batchmin; batch2 < batchmax; batch2++)
    {
      batch = 1 << batch2;
      timefft<cufftComplex>(n,batch,true, option);
    }
    
  for (batch2 = batchmin; batch2 < batchmax; batch2++)
    {
      batch = 1 << batch2;
      timefft<half2>(n,batch,false, option);
    }
  
  for (batch2 = batchmin; batch2 < batchmax; batch2++)
    {
      batch = 1 << batch2;
      timefft<half2>(n,batch,true, option);
    }
  
  fprintf(stdout, "Benchmark finished\n");

  return EXIT_SUCCESS;
}
