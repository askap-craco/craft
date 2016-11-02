#CUDA_PATH?=/Developer/NVIDIA/CUDA-7.5
CUDA_PATH=/usr/local/cuda-5.5/
objects=CudaTimer.o SigprocFileSet.o array.o boxcar.o cudafdmt.o fdmt.o fdmt_test.o fdmt_utils.o rescale.o CandidateSink.o CpuTimer.o DataSource.o SigprocFile.o
NVCC=$(CUDA_PATH)/bin/nvcc

NVCCFLAGS=-G -g -O3 -Xcompiler -fPIC --compile --relocatable-device-code=true -gencode arch=compute_20,code=compute_20 -gencode arch=compute_20,code=sm_20  -x cu 
NVLDFLAGS=--cudart static --relocatable-device-code=true -gencode arch=compute_20,code=compute_20 -gencode arch=compute_20,code=sm_20 -link

cudafdmt: $(objects)
	$(NVCC) $(NVLDFLAGS) -o cudafdmt $(objects)


$(objects): %o : %cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm cudafdmt $(objects)
