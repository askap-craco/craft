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
#define NBASELINE 276       // possible numbers: 276, 435, 435
//#define NBASELINE 1       // possible numbers: 276, 435, 435
#define NTIME_PER_BUFBLOCK 256    // configurable number, depends the available memory
#define NSAMP_PER_TIME (NCHAN*NBASELINE)

// Use float to emulate complex, the first float is real part and the second float is the image part
int prepare(float *in, float *cal, float *sky, float *out, float *average);
