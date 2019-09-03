/*******************************************************************************
** HOST Code HEADER FILE
*******************************************************************************/

#pragma once

#define CL_HPP_CL_1_2_DEFAULT_BUILD
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY 1
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS

#include <fcntl.h>
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdbool.h>
#include <CL/opencl.h>
#include <CL/cl_ext.h>

#define MAX_LENGTH 4096
#define MEM_ALIGNMENT 4096

#define CHANNELS  288       // possible numbers: 288, 288, 672
#define BASELINEs 276       // possible numbers: 276, 435, 435
#define SAMPLES_PER_BUF_BLOCK 256  // configurable number, depends the available memory
#define BUF_BLOCKS_PER_MODEL  8192 // Possible numbers: 32, 32, 64
#define MODELS 2            // number of models, defines the length of the try
