/*
 * cpu_kernels.h
 *
 *  Created on: 26 Sep 2016
 *      Author: ban115
 */

#ifndef CPU_KERNELS_H_
#define CPU_KERNELS_H_

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


#endif /* CPU_KERNELS_H_ */
