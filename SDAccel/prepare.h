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

#define NCHAN 288        // possible numbers: 288, 288, 672
#define NANT 30
//#define NBASELINE 2
#define NBASELINE (NANT*(NANT-1)/2)
//#define NTIME_PER_BUFBLOCK 2    // configurable number, depends the available memory
#define NTIME_PER_BUFBLOCK 64    // configurable number, depends the available memory
#define NSAMP_PER_TIME (NCHAN*NBASELINE)

// Use float to emulate complex, the first float is real part and the second float is the image part
int prepare(float *in, float *cal, float *sky, float *out, float *average);
