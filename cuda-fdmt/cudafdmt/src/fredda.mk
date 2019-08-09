# site-specific includes
# path to cuda - shoudl contain bin/nvcc
CUDA_PATH=/usr/

# target cuda architecture e.g. compute_60
CUDA_ARCH=compute_70

# target cuda code = e.g. sm_60
CUDA_CODE=sm_70


# path to cub.h from cub-1.6.4 https://nvlabs.github.io/cub/
CUB_PATH=$(KEITH_HOME)/build/cub/cub-1.6.4/cub/

# path to PSRDADA: http://psrdada.sourceforge.net/
#DADA_PATH=/home/craftop/askap/trunk/3rdParty/psrdada/psrdada-537159/install/
DADA_PATH=$(KEITH_HOME)/psrdada-install/
