.PHONY: help

help::
	$(ECHO) "Makefile Usage:"
	$(ECHO) "  make all TARGET=<sw_emu/hw_emu/hw> DEVICE=<FPGA platform> "
	$(ECHO) "      Command to pick the specific files and generates the design for specified Target and Device."
	$(ECHO) ""
	$(ECHO) "  make clean "
	$(ECHO) "      Command to remove the generated non-hardware files."
	$(ECHO) ""
	$(ECHO) "  make cleanall"
	$(ECHO) "      Command to remove all the generated files."
	$(ECHO) ""
	$(ECHO) "  make check TARGET=<sw_emu/hw_emu/hw> DEVICE=<FPGA platform>"
	$(ECHO) "      Command to run application in emulation."
	$(ECHO) ""

TARGET := sw_emu
DEVICE := xilinx_u280_xdma_201910_1
XCLBIN := ./xclbin
KERNEL_VADD := vadd

DSA := $(call device2sandsa, $(DEVICE))
BUILD_DIR := ./_x.$(TARGET)
BUILD_DIR_vadd = $(BUILD_DIR)/$(KERNEL_VADD)

SRC_DIR := .

# The C++ Compiler to use: XCPP 
CXX := $(XILINX_SDX)/bin/xcpp

# The kernel Compiler to use : XOCC
XOCC := $(XILINX_SDX)/bin/xocc

# The below are compile flags are passed to the C++ Compiler
opencl_CXXFLAGS += -g -I./ -I$(XILINX_XRT)/include -I$(XILINX_VIVADO)/include -Wall -O0 -g -std=c++14
# The below are linking flags for C++ Comnpiler 
opencl_LDFLAGS += -L$(XILINX_XRT)/lib -lOpenCL -lpthread

CXXFLAGS += $(opencl_CXXFLAGS)
LDFLAGS += $(opencl_LDFLAGS)

OBJECTS=util_sdaccel.o prepare.o

# VADD host CPP FILE
VADD_SRCS += $(SRC_DIR)/vadd_host.c $(SRC_DIR)/util_sdaccel.c
# VADD host Header FILE		    
VADD_HDRS += $(SRC_DIR)/vadd_host.h $(SRC_DIR)/util_sdaccel.h
# VADD executable name
VADD = vadd_host
EXECUTABLE += $(VADD)

# Other executable, which are for test of kernels
PREPARE_SRCS += $(SRC_DIR)/prepare_host.c $(SRC_DIR)/util_sdaccel.c $(SRC_DIR)/prepare.c
PREPARE_HDRS += $(SRC_DIR)/prepare_host.h $(SRC_DIR)/util_sdaccel.h $(SRC_DIR)/prepare.h
PREPARE = prepare_host
EXECUTABLE += $(PREPARE)

# Host compiler global settings
CXXFLAGS += -fmessage-length=0
LDFLAGS += -lrt -lstdc++ 

# Enable Profiling
REPORT := no
PROFILE:= yes
#timeline_trace is always enabled

# Kernel compiler global settings
CLFLAGS += -t $(TARGET) --platform $(DEVICE) --save-temps 

EMCONFIG_DIR = $(XCLBIN)

BINARY_CONTAINERS += $(XCLBIN)/$(KERNEL_VADD).$(TARGET).$(DEVICE).xclbin
BINARY_CONTAINER_vadd_OBJS += $(XCLBIN)/$(KERNEL_VADD).$(TARGET).$(DEVICE).xo

CP = cp -rf

.PHONY: exe
exe: $(EXECUTABLE)

.PHONY: all clean cleanall docs emconfig
all: check-devices exe $(BINARY_CONTAINERS) emconfig Makefile

# Building kernel
$(XCLBIN)/$(KERNEL_VADD).$(TARGET).$(DEVICE).xo: $(SRC_DIR)/vadd.cpp
	mkdir -p $(XCLBIN)
	$(XOCC) $(CLFLAGS) --temp_dir $(BUILD_DIR_vadd) -c -g -k $(KERNEL_VADD) -I'$(<D)' -o'$@' '$<'
$(XCLBIN)/$(KERNEL_VADD).$(TARGET).$(DEVICE).xclbin: $(BINARY_CONTAINER_vadd_OBJS)
	mkdir -p $(XCLBIN)
	$(XOCC) $(CLFLAGS) --temp_dir $(BUILD_DIR_vadd) -l $(LDCLFLAGS) --nk $(KERNEL_VADD):1 -o'$@' $(+)

$(VADD): $(VADD).o $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(VADD).o $(OBJECTS) -o $(VADD)

$(PREPARE): $(PREPARE).o $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(PREPARE).o $(OBJECTS) -o $(PREPARE)

%.o: %.c
	$(CXX) $(CXXFLAGS) -c $<

emconfig:$(EMCONFIG_DIR)/emconfig.json
$(EMCONFIG_DIR)/emconfig.json:
	emconfigutil --platform $(DEVICE) --od $(EMCONFIG_DIR)

check: all
ifeq ($(TARGET),$(filter $(TARGET),sw_emu hw_emu))
	$(CP) $(EMCONFIG_DIR)/emconfig.json .
	XCL_EMULATION_MODE=$(TARGET) ./$(VADD) $(XCLBIN)/$(KERNEL_VADD).$(TARGET).$(DEVICE).xclbin $(DEVICE)
else
	 ./$(VADD) $(XCLBIN)/$(KERNEL_VADD).$(TARGET).$(DEVICE).xclbin $(DEVICE)
endif
#	sdx_analyze profile -i profile_summary.csv -f html
#	sdx_analyze trace timeline_trace.csv
#	sdx -workspace workspace -report *.wdb
#


# Cleaning stuff
RMDIR = rm -rf 

clean:
	-$(RMDIR) $(EXECUTABLE) $(XCLBIN)/{*sw_emu*,*hw_emu*} 
	-$(RMDIR) sdaccel_* TempConfig system_estimate.xtxt *.rpt *.o *.csv
	-$(RMDIR) *.protoinst _xocc_* .Xil emconfig.json dltmp* xmltmp* *.log *.jou *.wcfg *.wdb

cleanall: clean
	-$(RMDIR) $(XCLBIN)
	-$(RMDIR) _x.*


check-devices:
ifndef DEVICE
	$(error DEVICE not set. Please set the DEVICE properly and rerun. Run "make help" for more details.)
endif

check_xrt:
ifndef XILINX_XRT
	$(error XILINX_XRT variable is not set, please set correctly and rerun)
endif


ECHO := @echo

#'estimate' for estimate report generation
#'system' for system report generation
ifneq ($(REPORT), no)
CLFLAGS += --report estimate
CLLDFLAGS += --report system
endif

#Generates profile summary report
ifeq ($(PROFILE), yes)
LDCLFLAGS += --profile_kernel data:all:all:all:all
LDCFLAGS += --profile_kernel  stall:all:all:all:all
LDCFALGS += --profile_kernel exec:all:all:all:all
endif
