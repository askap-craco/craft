#CUDA_PATH?=/Developer/NVIDIA/CUDA-7.5
#CUDA_PATH=/usr/local/cuda-5.5/

objects=CudaTimer.o SigprocFileSet.o array.o boxcar.o cudafdmt.o fdmt.o fdmt_utils.o rescale.o CandidateSink.o CpuTimer.o DataSource.o SigprocFile.o CandidateList.o DadaSource.o  InvalidSourceFormat.o Rescaler.o DataOrder.o DadaSet.o FilDirSet.o DadaSink.o Array4dDumper.o FreddaParams.o

fftbench_objects=CudaTimer.o fftbench.o CpuTimer.o
craftbench_objects=CudaTimer.o craftbench.o CpuTimer.o

all_objects=$(objects) fdmt_test.o fftbench.o craftbench.o 

include fredda.mk

NVCC=$(CUDA_PATH)/bin/nvcc

#Embed version
# See http://stackoverflow.com/questions/1704907/how-can-i-get-my-c-code-to-automatically-print-out-its-git-version-hash
GIT_VERSION := $(shell git describe --dirty --always --tags)

NVCCFLAGS=-g -O3  -Xcompiler -fopenmp -Xcompiler -fPIC -Xcompiler --compile --relocatable-device-code=true -gencode arch=$(CUDA_ARCH),code=$(CUDA_CODE)  -x cu -DVERSION=\"$(GIT_VERSION)\" -I$(CUB_PATH) -I$(DADA_PATH)/include 

NVLDFLAGS=--cudart static --relocatable-device-code=true -gencode arch=$(CUDA_ARCH),code=$(CUDA_CODE) -link -Xcompiler -L$(DADA_PATH)/lib -lpsrdada -lgomp -lcufft -g 


# Force rebuild on git change:http://stackoverflow.com/questions/3236145/force-gnu-make-to-rebuild-objects-affected-by-compiler-definition

all: git_version cudafdmt fftbenchapp fftbenchapp_deng

.PHONY: all

git_version: 
	echo '$(GIT_VERSION)' | cmp -s - $@ || echo '$(GIT_VERSION)' > $@

cudafdmt: $(objects) fdmt_test.o git_version
	$(NVCC) $(NVLDFLAGS) -o cudafdmt $(objects) fdmt_test.o 

fftbenchapp:  cufft_utils.h $(fftbench_objects)
	$(NVCC) $(NVLDFLAGS) -o fftbench $(fftbench_objects) 

fftbenchapp_deng:  cufft_utils.h $(craftbench_objects)
	$(NVCC) $(NVLDFLAGS) -o craftbench $(craftbench_objects) 

$(all_objects): %o  : %cu
	$(NVCC) $(NVCCFLAGS) -o $@ -c $<


.PHONY : clean

clean :
	rm cudafdmt fftbench $(objects)
