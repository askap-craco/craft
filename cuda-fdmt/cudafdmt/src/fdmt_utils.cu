/*
 * fdmtutils.cpp
 *
 *  Created on: 26 Oct 2016
 *      Author: ban115
 */

#include "fdmt_utils.h"

int __device__ calculate_offset()
{
	int ibeam = blockIdx.x;
	int nbeams = gridDim.x;

	int idt = blockIdx.y;
	int max_dt = gridDim.y;

	int nt = blockDim.x;
	int off = max_dt*(idt + ibeam*nbeams);
	return off;
}
