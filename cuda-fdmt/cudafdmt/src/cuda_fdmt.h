/*
 * cuda_fdmt.h
 *
 *  Created on: 4 Oct 2016
 *      Author: ban115
 */

#ifndef CUDA_FDMT_H_
#define CUDA_FDMT_H_

#include "fdmt.h"
#include "array.h"

void __global__ cuda_fdmt_iteration(fdmt_t fdmt,
                   int iteration_num,
                   array4d_t indata,
                   array4d_t outdata);




#endif /* CUDA_FDMT_H_ */
