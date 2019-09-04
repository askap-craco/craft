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
#define NTIME_PER_BUFBLOCK 256    // configurable number, depends the available memory
#define NTIME_AVERAGE 8192 // Possible numbers: 8192, 8192, 16384
#define NBUFBLOCK_AVERAGE (NTIME_AVERAGE/NTIME_PER_BUFBLOCK)
#define NAVEAGE 2                   // number of averages, defines the length of the try

// Use float to emulate complex, the first float is real part and the second float is the image part
int prepare(float *in, float *cal, float *sky, float *out, float *average);
