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

template <class intype>
int timefft(int n, int batch, bool inplace, int option)
{
  CudaTimer t;
  intype *data, *out_data, *data_host, *out_data_host;
  cufftHandle plan;
  size_t data_size=sizeof(intype)*n*(n/2 + 1)*batch;
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
  gpuErrchk(cudaMalloc((void**) &data, data_size));
  if (inplace)
    out_data = data;
  else
    gpuErrchk(cudaMalloc((void**) &out_data, data_size));

  /* Get memory on host */
  if(option != 0)
    {
      gpuErrchk(cudaMallocHost((void**) &data_host, data_size));
      if (inplace)
	out_data_host = data_host;
      else
	gpuErrchk(cudaMallocHost((void**) &out_data_host, data_size));
    }
  long long int nsize[] = {n,n };

  /* Create the FFT plan for multiple data types */
  size_t worksize;
  cufftSafeCall(cufftCreate(&plan));
  cufftSafeCall(cufftXtMakePlanMany(plan, 2, nsize,
				    NULL, 1, 0, itype,
				    NULL, 1, 0, otype,
				    batch, &worksize, etype
				    ));
  
  /* Warm up and do the real thing */
  cufftSafeCall(cufftXtExec(plan, data, out_data, CUFFT_INVERSE));
  int i, niter = 1000;
  for (i = 0; i < niter; ++i)
    {
      t.start();
      //cufftSafeCall(cufftExecC2R(plan, data, (outtype*) data));
      if(option == 1)
	gpuErrchk(cudaMemcpy(data, data_host, data_size, cudaMemcpyHostToDevice));
      cufftSafeCall(cufftXtExec(plan, data, out_data, CUFFT_INVERSE));
      if(option == 1)
	gpuErrchk(cudaMemcpy(out_data_host, out_data, data_size, cudaMemcpyDeviceToHost));
      t.stop(); 
    }

  /* Check the timer and display */
  float tavg_us = t.get_average_time() / float(batch) * 1e3f;
  fprintf(stdout, "%dx%d FFT batch=%d data=%d MB in-place=%d type=%d-> %d. Worksize=%d MB: %f microseconds/FFT= %f k FFTs/sec total=%0.2fs\n",
	  n,n,batch,data_size/1024/1024, inplace, itype,otype, worksize/1024/1024, tavg_us, 1./tavg_us*1e6f/1e3f);

  /* Destroy FFT plan */
  cufftSafeCall(cufftDestroy(plan));

  /* Free device memory */
  gpuErrchk(cudaFree(data));
  if (! inplace)
    gpuErrchk(cudaFree(out_data));

  /* Free host memory */
  if(option != 0)
    {
      gpuErrchk(cudaFreeHost(data_host));
      if (! inplace)
	gpuErrchk(cudaFreeHost(out_data_host));
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
