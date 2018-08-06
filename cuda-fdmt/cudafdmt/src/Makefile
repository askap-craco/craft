#CUDA_PATH?=/Developer/NVIDIA/CUDA-7.5
#CUDA_PATH=/usr/local/cuda-5.5/
CUDA_PATH=$(CUDA)
objects=CudaTimer.o SigprocFileSet.o array.o boxcar.o cudafdmt.o fdmt.o fdmt_test.o fdmt_utils.o rescale.o CandidateSink.o CpuTimer.o DataSource.o SigprocFile.o CandidateList.o
NVCC=$(CUDA_PATH)/bin/nvcc
CUB_PATH=$(HOME)/build/cub/cub-1.6.4/cub/

#Embed version
# See http://stackoverflow.com/questions/1704907/how-can-i-get-my-c-code-to-automatically-print-out-its-git-version-hash
GIT_VERSION := $(shell git describe --dirty --always --tags)

NVCCFLAGS=-g -O3 -Xcompiler -fPIC -Xcompiler -fopenmp --compile --relocatable-device-code=true -gencode arch=compute_60,code=sm_60  -x cu -DVERSION=\"$(GIT_VERSION)\" -I$(CUB_PATH)

NVLDFLAGS=--cudart static --relocatable-device-code=true -gencode arch=compute_60,code=compute_60 -gencode arch=compute_60,code=sm_60 -link -Xcompiler -fopenmp

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
