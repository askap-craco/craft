
COMMON = -Wall -g

CC=gcc $(COMMON)
CXX = g++ $(COMMON) -std=c++11 

%.o: %.cpp
	$(CXX) $(CCFLAGS) $< -c

%.o: %.cc
	$(CC) $(CXXFLAGS) $< -c

OBJ = CRAFTMode.o CodecCODIF.o CodecVCRAFT.o CodecVDIF1.o CommandLineParser.o DFH.o FileDescriptor.o FileProcessor.o SampleData.o VCRAFTBitMode.o VCRAFTParser.o main.o codifio.o dateutils.o

INCS = CRAFTMode.h CodecVDIF1.h FileProcessor.h VCRAFTParser.h		\
	Codec.h CommandLineParser.h SampleData.h codifio.h CodecCODIF.h DFH.h	\
	StdApp.h dateutils.h CodecVCRAFT.h FileDescriptor.h VCRAFTBitMode.h


all: $(OBJ) CRAFTConverter

CRAFTConverter: $(OBJ) $(INCS)
	$(CXX) -o CRAFTConverter $(OBJ)

clean:
	rm -rf *.o *~ $(OBJ) CRAFTConverter

.PHONY: clean all
