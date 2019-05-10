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
  CudaTimer tgpu;
  CpuTimer t;
  
	//intype *data, *out_data;
	cufftHandle plan;
	size_t data_size=sizeof(intype)*n*(n/2 + 1)*batch;
	int ngpus = 4;
	int gpus[]  = { 0,1,2,3 };
	/*
	for (int i = 0; i < ngpus; i++) {
	  gpuErrchk(cudaSetDevice(gpus[i]))
	  gpuErrchk(cudaMalloc((void**) &data, data_size));
	  if (inplace) {
	    out_data = data;
	  } else {
	    gpuErrchk(cudaMalloc((void**) &out_data, data_size));
	  }
	}
	*/

	cufftXtSubFormat format;
	if (inplace) {
	  format = CUFFT_XT_FORMAT_INPLACE;
	} else {
	  format = CUFFT_XT_FORMAT_INPUT; // THere is also _OUTPU
	}

	long long int nsize[] = {n,n };

		/*cufftSafeCall(cufftPlanMany(&plan, 2, n,
			NULL, 1, 0, // Simple input layout
			NULL, 1, 0, // Simple output layout
			CUFFT_C2R, BATCH));
	*/
	size_t worksize[4];
	cudaLibXtDesc* in_data;
	cudaLibXtDesc* out_data;
	cufftSafeCall(cufftCreate(&plan));
	if (ngpus > 1) {
	  cufftSafeCall(cufftXtSetGPUs(plan, ngpus, gpus));
	}
	///cufftSafeCall(cufftSetAutoAllocation());

	//cudaDataType itype = CUDA_C_16F;
	//cudaDataType etype = CUDA_C_16F;
	//cudaDataType otype = CUDA_R_16F;
	cufftSafeCall(cufftXtMakePlanMany(plan, 2, nsize,
					  NULL, 1, 256, itype,
					  NULL, 1, 256, otype,
					  batch, worksize, etype
					  ));

	cufftSafeCall(cufftXtMalloc(plan, (cudaLibXtDesc **)&in_data, format));
	if (inplace) {
	  out_data = in_data;
	} else {
	  cufftSafeCall(cufftXtMalloc(plan, (cudaLibXtDesc **)&out_data, format));
	}
	
	// warm up
	//cufftSafeCall(cufftXtExec(plan, in_data, out_data, CUFFT_INVERSE));
	cufftSafeCall(cufftXtExecDescriptorC2R(plan, in_data,  out_data));

	int niter = 100;
	t.start();
	tgpu.start();
	for (int i = 0; i < niter; ++i) {
	  cufftSafeCall(cufftXtExecDescriptorC2R(plan, in_data,  out_data));
		//cufftSafeCall(cufftExecC2R(plan, data, (outtype*) data));
		//cufftSafeCall(cufftXtExec(plan, in_data, out_data, CUFFT_INVERSE));
	}
	tgpu.stop();
	t.stop();
	cout << tgpu << endl;
	//cout << "CPU ONLY" << t << endl;

	float tavg_us = tgpu.get_average_time() / float(niter*batch) * 1e3f;

	printf("%dx%d FFT batch=%d data=%d MB in-place=%d type=%d-> %d. Worksize=%d MB: %f microseconds/FFT= %f k FFTs/sec\n",
	       n,n,batch,data_size/1024/1024, inplace, itype,otype, worksize[0]/1024/1024, tavg_us, 1./tavg_us*1e6f/1e3f);

	cufftSafeCall(cufftXtFree(in_data));
	if (! inplace) {
	  cufftSafeCall(cufftXtFree(out_data));
	}

	cufftSafeCall(cufftDestroy(plan));
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



