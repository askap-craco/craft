/*
 * cpu_kernels.h
 *
 *  Created on: 26 Sep 2016
 *      Author: ban115
 */

#ifndef CPU_KERNELS_H_
#define CPU_KERNELS_H_

#include "array.h"

void cpu_copy1(fdmt_dtype* dst, const fdmt_dtype* src, const int count) {
	assert(count >= 0);
	for(int i = 0; i < count; i++) {
		dst[i] = src[i];
	}
}

void cpu_copy2(fdmt_dtype* dst, const fdmt_dtype* src, const int count) {
	assert(count >= 0);
	memcpy(dst, src, count*sizeof(fdmt_dtype));
}

void cpu_sum1(fdmt_dtype* dst, const fdmt_dtype* src1, const fdmt_dtype* src2, const int count) {
	assert(count >= 0);
	for(int i = 0; i < count; i++) {
		dst[i] = src1[i] + src2[i];
	}
}

void array_cpu_copy1(array4d_t* dst, const array4d_t* src, coord3_t* dstidx, coord3_t* srcidx, int zcount)
{
	for(int w = 0; w < dst->nw; w++) {
		int srcidxi = array4d_idx(src, w, srcidx->x, srcidx->y, srcidx->z);
		int dstidxi = array4d_idx(dst, w, dstidx->x, dstidx->y, dstidx->z);
		for (int z = 0; z < zcount; z++) {
			dst->d[dstidxi + z] = src->d[srcidxi + z];
		}
	}
}

void array_cpu_sum1(array4d_t* dst, const array4d_t* src1, coord3_t* dstidx, coord3_t* src1idx, coord3_t* src2idx, int zcount)
{
	for(int w = 0; w < dst->nw; w++) {
		int src1idxi = array4d_idx(src1, w, src1idx->x, src1idx->y, src1idx->z);
		int src2idxi = array4d_idx(src1, w, src2idx->x, src2idx->y, src2idx->z);
		int dstidxi = array4d_idx(dst, w, dstidx->x, dstidx->y, dstidx->z);
		for (int z = 0; z < zcount; z++) {
			dst->d[dstidxi + z] = src1->d[src1idxi + z] + src1->d[src2idxi + z];
		}
	}
}


#endif /* CPU_KERNELS_H_ */
