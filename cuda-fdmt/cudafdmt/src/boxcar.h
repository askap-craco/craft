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
#include "CandidateList.h"

const int NBOX = 32; // Needs to be the warp size. One day we'll check that.
const int DT_BLOCKS =  1; // Number of DTs to do in a block simultaneously

int boxcar_do(array4d_t* indata, array4d_t* outdata);
int boxcar_do_cpu(const array4d_t* indata,
		array4d_t* outdata,
		array4d_t* boxcar_history,
		fdmt_dtype thresh, int max_ncand_per_block,  int mindm, int maxbc,
		CandidateSink& sink);

int boxcar_do_gpu(const array4d_t* indata,
		const array4d_t* weights,
		array4d_t* boxcar_data,
		array4d_t* boxcar_history,
		array4d_t* boxcar_discards,
		fdmt_dtype thresh, int max_ncand_per_block, int mindm, int maxbc,
		CandidateList* sink);

int boxcar_threshonly(const array4d_t* indata, size_t sampno, fdmt_dtype thresh, int max_ncand_per_block,int mindm,
		CandidateList& sink);

#endif /* BOXCAR_H_ */

