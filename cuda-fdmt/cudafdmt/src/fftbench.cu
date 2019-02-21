/*
 * fftbench.cu
 *
 *  Created on: 17 Sep 2018
 *      Author: ban115

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


//typedef half2 intype;
//typedef half outtype;

//typedef cufftComplex intype;
//typedef cufftReal outtype;
//typedef cufftComplex ftype;

template <class intype>
void timefft(int n, int batch, cudaDataType itype, cudaDataType etype, cudaDataType otype, bool inplace)
{
	CudaTimer t;
	intype *data, *out_data;
	cufftHandle plan;
	size_t data_size=sizeof(intype)*n*(n/2 + 1)*batch;
	gpuErrchk(cudaMalloc((void**) &data, data_size));
	if (inplace) {
	  out_data = data;
	} else {
	  gpuErrchk(cudaMalloc((void**) &out_data, data_size));
	}

	long long int nsize[] = {n,n };

		/*cufftSafeCall(cufftPlanMany(&plan, 2, n,
			NULL, 1, 0, // Simple input layout
			NULL, 1, 0, // Simple output layout
			CUFFT_C2R, BATCH));
	*/
	size_t worksize;
	cufftSafeCall(cufftCreate(&plan));
	//int gpus[]  = { cuda_device };
	//cufftSafeCall(cufftXtSetGPUs(plan, 1, gpus));
	//cufftSafeCall(cufftSetAutoAllocation());

	//cudaDataType itype = CUDA_C_16F;
	//cudaDataType etype = CUDA_C_16F;
	//cudaDataType otype = CUDA_R_16F;
	cufftSafeCall(cufftXtMakePlanMany(plan, 2, nsize,
			NULL, 1, 0, itype,
			NULL, 1, 0, otype,
			batch, &worksize, etype
			));

	// warm up
	cufftSafeCall(cufftXtExec(plan, data, out_data, CUFFT_INVERSE));

	int niter = 100;
	for (int i = 0; i < niter; ++i) {
		t.start();
		//cufftSafeCall(cufftExecC2R(plan, data, (outtype*) data));
		cufftSafeCall(cufftXtExec(plan, data, out_data, CUFFT_INVERSE));
		t.stop();

	}
	float tavg_us = t.get_average_time() / float(batch) * 1e3f;

	printf("%dx%d FFT batch=%d data=%d MB in-place=%d type=%d-> %d. Worksize=%d MB: %f microseconds/FFT= %f k FFTs/sec total=%0.2fs\n",
	       n,n,batch,data_size/1024/1024, inplace, itype,otype, worksize/1024/1024, tavg_us, 1./tavg_us*1e6f/1e3f);
	cufftSafeCall(cufftDestroy(plan));
	gpuErrchk(cudaFree(data));
	if (! inplace) {
	  gpuErrchk(cudaFree(out_data));
	}
}

int main(int argc, char* argv[])
{

	if (argc != 5) {
		printf("%s Usage: gpuid N batchmin batchmax\n", argv[0]);
		return EXIT_FAILURE;
	}

	int cuda_device = atoi(argv[1]);
	cudaDeviceProp p;
	gpuErrchk(cudaGetDeviceProperties(&p, cuda_device));

	printf("FFT Benchmark \n");
	printf("Device[%d]=%s v%d.%d Mem=%d GB shmem/block=%d constmem=%d Warp=%d Clock=%d MHz %d multiprocessors\n",
			cuda_device, p.name, p.major, p.minor, p.totalGlobalMem/1024/1024/1024, p.sharedMemPerBlock,
			p.totalConstMem, p.warpSize, p.clockRate/1000, p.multiProcessorCount
			);

	gpuErrchk( cudaSetDevice(cuda_device));

	int n = atoi(argv[2]);
	int batchmin = atoi(argv[3]);
	int batchmax = atoi(argv[4]);
	cudaDataType itype = CUDA_C_32F;
	cudaDataType etype = CUDA_C_32F;
	cudaDataType otype = CUDA_R_32F;

	for (int batch2 = batchmin; batch2 < batchmax; batch2++) {
		int batch = 1 << batch2;
		timefft<cufftComplex>(n,batch,itype, etype, otype, false);
	}

	for (int batch2 = batchmin; batch2 < batchmax; batch2++) {
		int batch = 1 << batch2;
		timefft<cufftComplex>(n,batch,itype, etype, otype, true);
	}


	itype = CUDA_C_16F;
	etype = CUDA_C_16F;
	otype = CUDA_R_16F;

	for (int batch2 = batchmin; batch2 < batchmax; batch2++) {
		int batch = 1 << batch2;
		timefft<half2>(n,batch,itype, etype, otype, false);
	}

	for (int batch2 = batchmin; batch2 < batchmax; batch2++) {
		int batch = 1 << batch2;
		timefft<half2>(n,batch,itype, etype, otype, true);
	}


	printf("Benchmark finished\n");
}



