/*******************************************************************************
** Kernel C Code HEADER FILE
*******************************************************************************/
#pragma once

#include <fcntl.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <complex>

typedef std::complex<float> cl_complex;
#define BURST_DATA_SIZE 8 // Pitch data into burst width, depends on the type of complex
#define BURST_DATA_HALF 4

#define NANT 4
//#define NCHAN 288        // possible numbers: 288, 288, 672, this minimum configuration is for the 512 bits width
#define NCHAN 256
#define NTIME_PER_BUFBLOCK 2    // configurable number, depends the available memory

//#define NANT 30
//#define NCHAN 288        // possible numbers: 288, 288, 672
////#define NCHAN 256        // possible numbers: 288, 288, 672
//#define NTIME_PER_BUFBLOCK 64    // configurable number, depends the available memory

#define NBASELINE (NANT*(NANT-1)/2)
#define NSAMP_PER_TIME (NCHAN*NBASELINE)

int prepare(cl_complex *in, cl_complex *cal, cl_complex *sky, cl_complex *out, cl_complex *average);

typedef struct burst_datatype {
  cl_complex data[BURST_DATA_SIZE];
} burst_datatype;

#define NBURST (NCHAN/BURST_DATA_SIZE)
