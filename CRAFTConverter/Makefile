
# If you don't have pkg-config setup uncomment the following and change as appropriate
#CODIFIO_ROOT=/set/me
#CODIFIO_INC=-I$(CODIFIO_ROOT)/include/codifio.h
#CODIFIO_LIB=-L$(CODIFIO_ROOT)/lib -lcodifio
CODIFIO_INC=$(shell pkg-config --cflags codifio)
CODIFIO_LIB=$(shell pkg-config --libs codifio)

COMMON = -Wall -Wunused-parameter -Wunused-variable -Wtype-limits -Wsign-compare -g $(CODIFIO_INC)

CC=gcc $(COMMON)
CXX = g++ $(COMMON) -std=c++11 

%.o: %.cpp
	$(CXX) $(CCFLAGS) $< -c

%.o: %.cc
	$(CC) $(CXXFLAGS) $< -c

OBJ = CRAFTMode.o CodecCODIF.o CodecVCRAFT.o CodecVDIF1.o CommandLineParser.o DFH.o FileDescriptor.o FileProcessor.o SampleData.o VCRAFTBitMode.o VCRAFTParser.o main.o dateutils.o

INCS = CRAFTMode.h CodecVDIF1.h FileProcessor.h VCRAFTParser.h		\
	Codec.h CommandLineParser.h SampleData.h CodecCODIF.h DFH.h	\
	StdApp.h dateutils.h CodecVCRAFT.h FileDescriptor.h VCRAFTBitMode.h

all: $(OBJ) CRAFTConverter


CRAFTConverter: $(OBJ) $(INCS)
	$(CXX) -o CRAFTConverter $(CODIFIO_LIB) $(OBJ)

clean:
	rm -rf *.o *~ $(OBJ) CRAFTConverter

.PHONY: clean all
