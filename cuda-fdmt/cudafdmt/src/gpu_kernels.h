/*
 * gpu_kernels.h
 *
 *  Created on: 27 Sep 2016
 *      Author: ban115
 */

#ifndef GPU_KERNELS_H_
#define GPU_KERNELS_H_

#include <cuda.h>

// Simple delay & sum over multiple dimensions
// Assumes BFDT ordering
__global__ void delay_sum1(float* in, float* out, int dt)
{
	int beam = blockIdx.x;
	int nbeams = gridDim.x;
	int chan = blockIdx.y;
	int nchans = gridDim.y;
}

__global__ void _array_gpu_copy_kernel(fdmt_dtype* dst, fdmt_dtype* src, int dst_stride, int src_stride)
{
	int beam = blockIdx.x;
	int t = threadIdx.x;
	int srcidxi = beam*src_stride + t;
	int dstidxi = beam*dst_stride + t;
	dst[dstidxi] = src[srcidxi];
}

__global__ void _array_gpu_sum_kernel(fdmt_dtype* dst, fdmt_dtype* src1, fdmt_dtype* src2, int dst_stride, int src_stride)
{
	int beam = blockIdx.x;
	int t = threadIdx.x;
	int srcidxi = beam*src_stride + t;
	int dstidxi = beam*dst_stride + t;
	dst[dstidxi] = src1[srcidxi] + src2[srcidxi];
}

__host__ void array_gpu_copy1(array4d_t* dst, const array4d_t* src, coord3_t* dstidx, coord3_t* srcidx, int zcount)
{
	int src_offset = array4d_idx(src, 0, srcidx->x, srcidx->y, srcidx->z);
	int src_stride = array4d_idx(src, 1, 0, 0, 0);
	int dst_offset = array4d_idx(dst, 0, dstidx->x, dstidx->y, dstidx->z);
	int dst_stride = array4d_idx(dst, 1, 0, 0, 0);

	int nbeams = dst->nw;

	_array_gpu_copy_kernel<<<nbeams, zcount>>>(dst->d_device + dst_offset, src->d_device + src_offset, dst_stride, src_stride);
}

__host__  void array_gpu_sum1(array4d_t* dst, const array4d_t* src, coord3_t* dstidx, coord3_t* src1idx, coord3_t* src2idx, int zcount)
{
	int src1_offset = array4d_idx(src, 0, src1idx->x, src1idx->y, src1idx->z);
	int src2_offset = array4d_idx(src, 0, src2idx->x, src2idx->y, src2idx->z);
	int src_stride = array4d_idx(src, 1, 0, 0, 0);

	int dst_offset = array4d_idx(dst, 0, dstidx->x, dstidx->y, dstidx->z);
	int dst_stride = array4d_idx(dst, 1, 0, 0, 0);

	int nbeams = dst->nw;

	_array_gpu_sum_kernel<<<nbeams, zcount>>>(dst->d_device + dst_offset, src->d_device + src1_offset, src->d_device + src2_offset, dst_stride, src_stride);

}


#endif /* GPU_KERNELS_H_ */
