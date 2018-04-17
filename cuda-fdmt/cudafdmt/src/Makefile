#CUDA_PATH?=/Developer/NVIDIA/CUDA-7.5
#CUDA_PATH=/usr/local/cuda-5.5/

objects=CudaTimer.o SigprocFileSet.o array.o boxcar.o cudafdmt.o fdmt.o fdmt_test.o fdmt_utils.o rescale.o CandidateSink.o CpuTimer.o DataSource.o SigprocFile.o CandidateList.o DadaSource.o  InvalidSourceFormat.o Rescaler.o DataOrder.o


include fredda.mk

NVCC=$(CUDA_PATH)/bin/nvcc

#Embed version
# See http://stackoverflow.com/questions/1704907/how-can-i-get-my-c-code-to-automatically-print-out-its-git-version-hash
GIT_VERSION := $(shell git describe --dirty --always --tags)

NVCCFLAGS=-g -O0 -Xcompiler -fPIC -Xcompiler -fopenmp --compile --relocatable-device-code=true -gencode arch=$(CUDA_ARCH),code=$(CUDA_CODE)  -x cu -DVERSION=\"$(GIT_VERSION)\" -I$(CUB_PATH) -I$(DADA_PATH)/include -std=c++11

NVLDFLAGS=--cudart static --relocatable-device-code=true -gencode arch=$(CUDA_ARCH),code=$(CUDA_CODE) -gencode arch=$(CUDA_ARCH),code=$(CUDA_CODE) -link -Xcompiler -fopenmp -L$(DADA_PATH)/lib -lpsrdada -g


# Force rebuild on git change:http://stackoverflow.com/questions/3236145/force-gnu-make-to-rebuild-objects-affected-by-compiler-definition

all: git_version cudafdmt

.PHONY: all

git_version: 
	echo '$(GIT_VERSION)' | cmp -s - $@ || echo '$(GIT_VERSION)' > $@

cudafdmt: $(objects) git_version
	$(NVCC) $(NVLDFLAGS) -o cudafdmt $(objects)

$(objects): %o  : %cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm cudafdmt $(objects)
