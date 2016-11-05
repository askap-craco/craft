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

int boxcar_do(array4d_t* indata, array4d_t* outdata);
int boxcar_threshonly(const array4d_t* indata, fdmt_dtype thresh,
		CandidateSink& sink);
#endif /* BOXCAR_H_ */

