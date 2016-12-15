/*
 * boxcar.h
 *
 *  Created on: 26 Oct 2016
 *      Author: ban115
 */

#ifndef BOXCAR_H_
#define BOXCAR_H_

#include "array.h"
#include "CandidateSink.h"

const int NBOX = 32; // Needs to be the warp size. One day we'll check that.

int boxcar_do(array4d_t* indata, array4d_t* outdata);
int boxcar_do_cpu(const array4d_t* indata,
		array4d_t* outdata,
		array4d_t* boxcar_history,
		size_t sampno,
		fdmt_dtype thresh, int max_ncand_per_block, int mindm,
		CandidateSink& sink);
int boxcar_threshonly(const array4d_t* indata, size_t sampno, fdmt_dtype thresh, int max_ncand_per_block,int mindm,
		CandidateSink& sink);
#endif /* BOXCAR_H_ */

