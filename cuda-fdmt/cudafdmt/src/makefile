DEBUG ?= 1
ifeq ($(DEBUG), 1)
    C_FLAGS = -DDEBUG
else
    C_FLAGS = -DNDEBUG
endif

auxiliary_objects = CudaTimer.o CpuTimer.o
main_objects      = fftbench.o fftbench_1d.o fftbench_callback.o
all_objects       = $(auxiliary_objects) $(main_objects)

all: fftbench fftbench_1d fftbench_callback

C_FLAGS      += -g
CU_FLAGS     = -dc -ccbin g++ -m64  \
	-gencode=arch=compute_30,code=sm_30 \
	-gencode=arch=compute_50,code=sm_50 \
	-gencode=arch=compute_52,code=sm_52 \
	-gencode=arch=compute_60,code=sm_60 \
	-gencode=arch=compute_61,code=sm_61 \
	-gencode=arch=compute_70,code=sm_70 \
	-gencode=arch=compute_70,code=compute_70 \
	-O2

NVCC         = nvcc
SRC_DIR      = .
OBJ_DIR      = .

LIBS         = -lcudart -lcuda -lcufft
LIBS_STATIC  = -lcudart -lcuda -lcufft_static -lculibos
LIB_DIRS     = 
INCLUDE_DIRS = 

fftbench_callback:fftbench_callback.o $(auxiliary_objects)
	$(NVCC) -ccbin g++ -m64 -o fftbench_callback fftbench_callback.o $(auxiliary_objects) $(LIB_DIRS) $(LIBS_STATIC)

fftbench:fftbench.o $(auxiliary_objects)
	$(NVCC) -ccbin g++ -m64 -o fftbench fftbench.o $(auxiliary_objects) $(LIB_DIRS) $(LIBS_STATIC)

fftbench_1d:fftbench_1d.o $(auxiliary_objects)
	$(NVCC) -ccbin g++ -m64 -o fftbench_1d fftbench_1d.o $(auxiliary_objects) $(LIB_DIRS) $(LIBS)

$(all_objects): %o  : %cu
	$(NVCC) $(INCLUDE_DIRS) ${C_FLAGS} $(CU_FLAGS) -o $@ -c $<

clean:
	rm -f *.o *~
