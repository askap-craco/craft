DEBUG ?= 1
ifeq ($(DEBUG), 1)
    C_FLAGS = -DDEBUG
else
    C_FLAGS = -DNDEBUG
endif

C_FLAGS      += -g
#CU_FLAGS     = -rdc=true -Wno-deprecated-gpu-targets -gencode arch=compute_20,code=sm_20 -gencode arch=compute_20,code=sm_21 -gencode arch=compute_30,code=sm_30 -gencode arch=compute_35,code=sm_35 -gencode arch=compute_50,code=sm_50 -gencode arch=compute_52,code=sm_52
#CU_FLAGS     = -rdc=true -Wno-deprecated-gpu-targets --default-stream per-thread
#CU_FLAGS     = -Wno-deprecated-gpu-targets --default-stream per-thread -dc -ccbin g++ -m64 -Xcompiler -static-libgcc -Xcompiler -static-libstdc++ -gencode=arch=compute_50,code=sm_50 \
#  -gencode=arch=compute_52,code=sm_52 \
#  -gencode=arch=compute_60,code=sm_60 \
#  -gencode=arch=compute_61,code=sm_61 \
#  -gencode=arch=compute_70,code=sm_70 \
#  -gencode=arch=compute_70,code=compute_70 \
#  -O2

CU_FLAGS     = -dc -ccbin g++ -m64   -gencode=arch=compute_30,code=sm_30 \
-gencode=arch=compute_50,code=sm_50 \
  -gencode=arch=compute_52,code=sm_52 \
  -gencode=arch=compute_60,code=sm_60 \
  -gencode=arch=compute_61,code=sm_61 \
  -gencode=arch=compute_70,code=sm_70 \
  -gencode=arch=compute_70,code=compute_70 \
  -O2

include fredda.mk

#NVCC         = $(CUDA_PATH)/bin/nvcc
NVCC         = nvcc
SRC_DIR      = .
OBJ_DIR      = .

LIBS         = -lpsrdada -lcudart -lcuda -lm -lrt -lcufft -lpthread 
LIBS_STATIC  = -lpsrdada -lcudart -lcuda -lm -lrt -lpthread -lcufft_static -lculibos
LIB_DIRS     = -L/usr/local/cuda/lib64 -L/usr/local/lib -L$(DADA_PATH)/lib -L$/usr/lib/x86_64-linux-gnu/ 
INCLUDE_DIRS = -I/usr/local/include -I$(CUB_PATH) -I$(DADA_PATH)/include -I/usr/include

all: fftbench fftbench_1d fftbench_callback 

fftbench_callback:fftbench_callback.o CudaTimer.o CpuTimer.o
	$(NVCC) -ccbin g++ -m64 -o fftbench_callback fftbench_callback.o CudaTimer.o CpuTimer.o $(LIB_DIRS) $(LIBS_STATIC)

fftbench_callback.o:fftbench_callback.cu
	$(NVCC) -c fftbench_callback.cu $(INCLUDE_DIRS) ${C_FLAGS} ${CU_FLAGS}

fftbench:fftbench.o CudaTimer.o CpuTimer.o
	$(NVCC) -ccbin g++ -m64 -o fftbench fftbench.o CudaTimer.o CpuTimer.o $(LIB_DIRS) $(LIBS)

fftbench.o:fftbench.cu
	$(NVCC) -c fftbench.cu $(INCLUDE_DIRS) ${C_FLAGS} ${CU_FLAGS}

fftbench_1d:fftbench_1d.o CudaTimer.o CpuTimer.o
	$(NVCC) -ccbin g++ -m64 -o fftbench_1d fftbench_1d.o CudaTimer.o CpuTimer.o $(LIB_DIRS) $(LIBS)

fftbench_1d.o:fftbench_1d.cu
	$(NVCC) -c fftbench_1d.cu $(INCLUDE_DIRS) ${C_FLAGS} ${CU_FLAGS}

CudaTimer.o:CudaTimer.cu
	$(NVCC) -c CudaTimer.cu $(INCLUDE_DIRS) ${C_FLAGS} ${CU_FLAGS}

CpuTimer.o:CpuTimer.cu
	$(NVCC) -c CpuTimer.cu $(INCLUDE_DIRS) ${C_FLAGS} ${CU_FLAGS}

clean:
	rm -f *.o *~
