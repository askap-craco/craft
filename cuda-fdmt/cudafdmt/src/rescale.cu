/*
 * Rescaling utilities
 * Author: Keith Bannister <keith.bannister@csiro.au>
 */

#include "rescale.h"
#include <stdint.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include "array.h"

typedef float rescale_dtype;

void* rescale_malloc(size_t sz)
{
	void* ptr = malloc(sz);
	assert(ptr);
	return ptr;
}

void rescale_update_scaleoffset(rescale_t* rescale)
{
	//assert(rescale->interval_samps >= 0);
	assert(rescale->target_stdev > 0);
	float nsamp = (float) rescale->sampnum;
	for (unsigned i = 0; i < rescale->num_elements; i++) {
		float mean = rescale->sum[i]/nsamp;
		float meansq = rescale->sumsq[i]/nsamp;
		float variance = meansq - mean*mean;

		if (rescale->interval_samps == 0) { // Don't do rescaling
			rescale->scale[i] = 1.0;
			rescale->offset[i] = 0.0;
		} else {

			if (variance == 0.0) {
				rescale->scale[i] = rescale->target_stdev;
			} else {
				rescale->scale[i] = rescale->target_stdev / sqrt(variance);
			}

			rescale->offset[i] = -mean + rescale->target_mean/rescale->scale[i];
		}

		// reset values to zero
		rescale->sum[i] = 0.0;
		rescale->sumsq[i] = 0.0;
	}

	rescale->sampnum = 0;
}

void rescale_update_none(rescale_t* rescale, float* inx, float*outx) 
{
	for (unsigned i = 0; i <  rescale->num_elements; i++) {
		float vin = inx[i];
		outx[i] = vin;
	}
}


void rescale_update_float(rescale_t* rescale, float* fdata, float* sampbuf, unsigned istart)
{
	float* inx = &fdata[istart];
	float* outx = &sampbuf[istart];
	for (unsigned i = 0; i <  rescale->num_elements; i++) {
		float vin = inx[i];
		float vin2 = vin*vin;
		rescale->sum[i] += vin;
		rescale->sumsq[i] += vin2;
		outx[i] = (vin + rescale->offset[i]) * rescale->scale[i];
	}
}

void rescale_update_float_polsum(rescale_t* rescale, float* fdata, float* sampbuf, unsigned istart)
{
	float* inx = &fdata[istart];
	float* outx = &sampbuf[istart];

	for (unsigned i = 0; i < rescale->num_elements/2; i++) {
		unsigned j = 2*i;
		float vin = inx[j];
		rescale->sum[j] += vin;
		rescale->sumsq[j] += vin*vin;

		float uin = inx[j+1];
		rescale->sum[j+1] += uin;
		rescale->sumsq[j+1] += uin*uin;

		float vscale = (vin + rescale->offset[j]) * rescale->scale[j];
		float uscale = (uin + rescale->offset[j+1]) * rescale->scale[j+1];

		float vout = (vscale + uscale)/2.0;

		outx[j] = vout;
	}
}

void rescale_update_uint8(rescale_t* rescale, float* fdata, uint8_t* sampbuf, unsigned istart)
{
	float* in = &fdata[istart];
	uint8_t* out = &sampbuf[istart];

	for (unsigned i = 0; i < rescale->num_elements; i++) {
		float vin = in[i];
		rescale->sum[i] += vin;
		rescale->sumsq[i] += vin*vin;
		float vout = (vin + rescale->offset[i]) * rescale->scale[i];

		if (vout < 0) {
			out[i] = 0;
		} else if (vout > 255) {
			out[i] = 255;
		} else {
			out[i] = (uint8_t) vout;
		}
	}
}

void rescale_update_uint8_polsum(rescale_t* rescale, float* fdata, uint8_t* sampbuf, unsigned istart)
{
	float* in = &fdata[istart];
	uint8_t* out = &sampbuf[istart];

	for (unsigned i = 0; i < rescale->num_elements/2; i++) {
		unsigned j=2*i;
		float vin = in[j];
		rescale->sum[j] += vin;
		rescale->sumsq[j] += vin*vin;

		float uin=in[j+1];
		rescale->sum[j+1] += uin;
		rescale->sumsq[j+1] += uin*uin;

		float vscale = (vin + rescale->offset[j]) * rescale->scale[j];
		float uscale = (uin + rescale->offset[j+1]) * rescale->scale[j+1];

		float vout = (vscale+uscale)/2.0;
		if (vout < 0) {
			out[j] = 0;
		} else if (vout > 255) {
			out[j] = 255;
		} else {
			out[j] = (uint8_t) vout;
		}
	}
}

void rescale_update_int8(rescale_t* rescale, float* __restrict__ in, int8_t* __restrict__ out)
{

	for (unsigned i = 0; i < rescale->num_elements; i++) {
		float vin = in[i];
		rescale->sum[i] += vin;
		rescale->sumsq[i] += vin*vin;
		float vout = (vin + rescale->offset[i]) * rescale->scale[i];
		if (vout < -128) {
			out[i] = -128;
		} else if (vout > 127) {
			out[i] = 127;
		} else {
			out[i] = (int8_t) vout;
		}
	}
	rescale->sampnum++;
	if (rescale->sampnum >= rescale->interval_samps) {
		rescale_update_scaleoffset(rescale);
	}

}

float rescale_update_decay_float_single(rescale_t* rescale, uint64_t i, float vin)
{
	assert(i < rescale->num_elements);
	float k = rescale->decay_constant;
	assert(k >= 0);

	rescale->sum[i] += vin;
	rescale->sumsq[i] += vin*vin;
	float vout = (vin + rescale->offset[i]) * rescale->scale[i];
	rescale->decay_offset[i] = (vout + rescale->decay_offset[i]*k) / (1.0 + k);

	float out = vout - rescale->decay_offset[i];

	return out;
}

void rescale_update_decay_float(rescale_t* rescale, float* __restrict__ in, float* __restrict__ out)
{
	float k = rescale->decay_constant;

	for (unsigned i = 0; i < rescale->num_elements; i++) {
		float vin = in[i];
		rescale->sum[i] += vin;
		rescale->sumsq[i] += vin*vin;
		float vout = (vin + rescale->offset[i]) * rescale->scale[i];
		rescale->decay_offset[i] = (vout + rescale->decay_offset[i]*k) / (1.0 + k);
		out[i] = vout - rescale->decay_offset[i];
	}
	rescale->sampnum++;
	if (rescale->sampnum >= rescale->interval_samps) {
		rescale_update_scaleoffset(rescale);
	}

}

void rescale_update_decay_uint8(rescale_t* rescale,  float* in, uint8_t* out)
{
	float k = rescale->decay_constant;

	for (unsigned i = 0; i < rescale->num_elements; i++) {
		float vin = in[i];
		rescale->sum[i] += vin;
		rescale->sumsq[i] += vin*vin;
		float vout = (vin + rescale->offset[i]) * rescale->scale[i];
		rescale->decay_offset[i] = (vout + rescale->decay_offset[i]*k) / (1.0 + k);
		float rout = (vout - rescale->decay_offset[i]);
		if (rout < 0) {
			out[i] = 0;
		} else if (rout > 255) {
			out[i] = 255;
		} else {
			out[i] = (uint8_t) rout;
		}

	}

	rescale->sampnum++;
	if (rescale->sampnum >= rescale->interval_samps) {
		rescale_update_scaleoffset(rescale);
	}
}

rescale_t* rescale_allocate(rescale_t* rescale, uint64_t nelements) 
{
	size_t sz = nelements*sizeof(float);

	rescale->sum = (float*) rescale_malloc(sz);
	rescale->sumsq = (float*) rescale_malloc(sz);
	rescale->scale = (float*) rescale_malloc(sz);
	rescale->offset = (float*) rescale_malloc(sz);
	rescale->decay_offset = (float*) rescale_malloc(sz);
	rescale->sampnum = 0;
	rescale->num_elements = nelements;

	for(uint64_t i = 0; i < nelements; ++i) {
		rescale->sum[i] = 0;
		rescale->sumsq[i] = 0;
		rescale->scale[i] = 1.0;
		rescale->offset[i] = 0.0;
		rescale->decay_offset[i] = 0.0;
	}

	//rescale_update_scaleoffset(rescale);
	return rescale;

}

rescale_dtype* rescale_cumalloc(uint64_t sz)
{
	rescale_dtype* ptr;
	gpuErrchk(cudaMalloc((void**) &ptr, sz));
	gpuErrchk(cudaMemset(ptr, 0, sz));
	return ptr;
}

void rescale_arraymalloc(array4d_t* arr, uint64_t sz)
{
	arr->nw = 1;
	arr->nx = 1;
	arr->ny = 1;
	arr->nz = sz;
	array4d_malloc(arr);
	array4d_set(arr, 0);
}

rescale_gpu_t* rescale_allocate_gpu(rescale_gpu_t* rescale, uint64_t nelements)
{
	size_t sz = nelements*sizeof(rescale_dtype);
	rescale_arraymalloc(&rescale->sum, nelements);
	rescale_arraymalloc(&rescale->sumsq, nelements);
	rescale_arraymalloc(&rescale->scale, nelements);
	rescale_arraymalloc(&rescale->offset, nelements);
	rescale_arraymalloc(&rescale->decay_offset, nelements);
	array4d_set(&rescale->scale, 1.0);

	rescale->sampnum = 0;
	rescale->num_elements = nelements;

	return rescale;

}

void rescale_set_scale_offset_gpu(rescale_gpu_t* rescale, float scale, float offset)
{
	array4d_set(&rescale->scale, scale);
	array4d_set(&rescale->offset, offset);
}



void rescale_update_and_transpose_float(rescale_t& rescale, array4d_t& read_arr, array4d_t& rescale_buf, uint8_t* read_buf, bool invert_freq)
{
	int nbeams = rescale_buf.nw;
	int nf = rescale_buf.nx;
	int nt = rescale_buf.nz;
	assert(rescale_buf.ny == 1);
	rescale.sampnum += nt;
	for(int t = 0; t < nt; ++t) {
#pragma omp parallel for
		for (int b = 0; b < nbeams; ++b) {
			int instart = array4d_idx(&read_arr, 0, b, t, 0);

			for (int f = 0; f < nf; ++f) {
				// NOTE: FDMT expects channel[0] at fmin
				// so invert the frequency axis if the frequency offset is negative
				int outf = f;
				if (invert_freq) {
					outf = nf - f - 1;
				}
				int inidx = instart + f;
				int outidx = array4d_idx(&rescale_buf, b, outf, 0, t);

				//printf("t=%d b=%d f=%d inidx=%d outidx=%d\n", t, b, f, inidx, outidx);
				// writes to inbuf
				size_t rs_idx = outf + nf*b;
				float v_rescale;
				//printf("Rescaling to mean=%f stdev=%f decay constant=%f\n",rescale.target_mean,rescale.target_stdev, rescale.decay_constant);

				v_rescale = rescale_update_decay_float_single(&rescale, rs_idx, (float) read_buf[inidx]);
				rescale_buf.d[outidx] = v_rescale;
				//printf("block=%d t=%d b=%d f=%d vin=%d vout=%f \n", blocknum, t, b, f, read_buf[inidx], v_rescale);

			}
		}
	}
}


__global__ void rescale_update_and_transpose_float_kernel(
		const uint8_t* __restrict__ inarr,
		rescale_dtype* __restrict__ sumarr,
		rescale_dtype* __restrict__ sumsqarr,
		rescale_dtype* __restrict__ decay_offsetarr,
		const rescale_dtype* __restrict__ offsetarr,
		const rescale_dtype* __restrict__ scalearr,
		rescale_dtype* __restrict__ outarr,
		float decay_constant,
		int nt,
		bool invert_freq)
{
	int ibeam = blockIdx.x;
	int c = threadIdx.x;
	int nf = blockDim.x;
	const rescale_dtype k = decay_constant;

	int rsidx = c + nf*ibeam;
	// all these reads are nice and coalesced
	rescale_dtype sum = sumarr[rsidx]; // read from global memory
	rescale_dtype sumsq = sumsqarr[rsidx]; // read from globa.
	rescale_dtype decay_offset = decay_offsetarr[rsidx];  // read from globa.
	rescale_dtype offset = offsetarr[rsidx]; // read from global
	rescale_dtype scale = scalearr[rsidx]; // read from global

	int outc;
	if (invert_freq) {
		outc = nf - 1 - c;
	} else {
		outc = c;
	}

	// input = BTF order
	// output = BFT order
	// Rescale order: BF
	for (int t = 0; t < nt; ++t) {
		int inidx = c + nf*(t + nt*ibeam);
		int outidx = t + nt*(outc + nf*ibeam);
		// coalesced read from global
		rescale_dtype vin = (rescale_dtype)inarr[inidx]; // read from global
		sum += vin;
		sumsq +=  vin*vin;
		rescale_dtype vout = (vin + offset) * scale;
		decay_offset = (vout + decay_offset*k)/(1.0 + k);

		// non-coalesced write (transpose. Sorry)
		rescale_dtype sout = vout - decay_offset;
		outarr[outidx] = sout;

	}

	// write everything back to global memory -- all coalesced
	sumarr[rsidx] = sum;
	sumsqarr[rsidx] = sumsq;
	decay_offsetarr[rsidx] = decay_offset;

}

void rescale_update_and_transpose_float_gpu(rescale_gpu_t& rescale, array4d_t& rescale_buf, const uint8_t* read_buf, bool invert_freq)
{
	int nbeams = rescale_buf.nw;
	int nf = rescale_buf.nx;
	int nt = rescale_buf.nz;
	assert(rescale_buf.ny == 1);
	rescale_update_and_transpose_float_kernel<<<nbeams, nf>>>(
			read_buf,
			rescale.sum.d_device,
			rescale.sumsq.d_device,
			rescale.decay_offset.d_device,
			rescale.offset.d_device,
			rescale.scale.d_device,
			rescale_buf.d_device,
			rescale.decay_constant,
			nt,
			invert_freq);
	rescale.sampnum += nt;
	gpuErrchk(cudaDeviceSynchronize());

}

__global__ void rescale_update_scaleoffset_kernel(
		rescale_dtype* __restrict__ sum,
		rescale_dtype* __restrict__ sumsq,
		rescale_dtype* __restrict__ offset,
		rescale_dtype* __restrict__ scalearr,
		rescale_dtype nsamp,
		rescale_dtype target_stdev,
		rescale_dtype target_mean)
{
	int c = threadIdx.x;
	int nf = blockDim.x;
	int ibeam = blockIdx.x;
	int i = c + nf*ibeam;
	rescale_dtype mean = sum[i]/nsamp;
	rescale_dtype meansq = sumsq[i]/nsamp;
	rescale_dtype variance = meansq - mean*mean;
	rescale_dtype scale;

	if (variance == 0.0) {
		scale = target_stdev;
	} else {
		scale = target_stdev / sqrt(variance);
	}

	offset[i] = -mean + target_mean/scale;
	scalearr[i] = scale;

	// reset values to zero
	sum[i] = 0.0;
	sumsq[i] = 0.0;
}
void rescale_update_scaleoffset_gpu(rescale_gpu_t& rescale)
{
	assert(rescale.interval_samps > 0);
	int nthreads = 336;
	assert(rescale.num_elements % nthreads == 0);
	int nblocks = rescale.num_elements / nthreads;
	rescale_update_scaleoffset_kernel<<<nblocks, nthreads>>>(
			rescale.sum.d_device,
			rescale.sumsq.d_device,
			rescale.offset.d_device,
			rescale.scale.d_device,
			(rescale_dtype) rescale.sampnum,
			rescale.target_stdev,
			rescale.target_mean);
	gpuErrchk(cudaDeviceSynchronize());
	rescale.sampnum = 0;
}

