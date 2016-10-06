CUDA_PATH?=/Developer/NVIDIA/CUDA-7.5
objects=CudaTimer.o array.o cuda_fdmt.o cudafdmt.o fdmt.o fdmt_test.o
NVCC=$(CUDA_PATH)/bin/nvcc
NVCCFLAGS=-G -g -O3 -Xcompiler -fPIC --compile --relocatable-device-code=true -gencode arch=compute_30,code=compute_30 -gencode arch=compute_30,code=sm_30  -x cu
NVLDFLAGS=--cudart static --relocatable-device-code=true -gencode arch=compute_30,code=compute_30 -gencode arch=compute_30,code=sm_30 -link -o

cudafdmt: $(objects)
	$(NVCC) $(NVLDFLAFGS) -o cudafdmt $(objects)


$(objects): %o : %cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm cudafdmt $(objects)
