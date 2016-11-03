#CUDA_PATH?=/Developer/NVIDIA/CUDA-7.5
CUDA_PATH=/usr/local/cuda-5.5/
objects=CudaTimer.o SigprocFileSet.o array.o boxcar.o cudafdmt.o fdmt.o fdmt_test.o fdmt_utils.o rescale.o CandidateSink.o CpuTimer.o DataSource.o SigprocFile.o
NVCC=$(CUDA_PATH)/bin/nvcc

#Embed version
# See http://stackoverflow.com/questions/1704907/how-can-i-get-my-c-code-to-automatically-print-out-its-git-version-hash
GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)

NVCCFLAGS=-G -g -O0 -Xcompiler -fPIC -Xcompiler -fopenmp --compile --relocatable-device-code=true -gencode arch=compute_20,code=compute_20 -gencode arch=compute_20,code=sm_20  -x cu -DVERSION=\"$(GIT_VERSION)\"

NVLDFLAGS=--cudart static --relocatable-device-code=true -gencode arch=compute_20,code=compute_20 -gencode arch=compute_20,code=sm_20 -link -Xcompiler -fopenmp

# Force rebuild on git change:http://stackoverflow.com/questions/3236145/force-gnu-make-to-rebuild-objects-affected-by-compiler-definition

#.PHONY: force
#git_version: force
#	echo '$(GIT_VERSION)' | cmp -s - $@ || echo '$(GIT_VERSION)' > $@

cudafdmt: $(objects)
	$(NVCC) $(NVLDFLAGS) -o cudafdmt $(objects)

$(objects): %o  : %cu
	$(NVCC) $(NVCCFLAGS) -o $@ $<

.PHONY : clean

clean :
	rm cudafdmt $(objects)
