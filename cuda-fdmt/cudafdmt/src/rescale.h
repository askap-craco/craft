#ifndef _RESCALE_H
#define _RESCALE_H

#include <stdint.h>
#include "array.h"

typedef float rescale_dtype;

typedef struct _rescale_t{
	/* stuff for rescaling */
	float* sum;
	float* sum2;
	float* scale;
	float* offset;
	float* decay_offset;
	uint64_t interval_samps;
	uint64_t sampnum;
	int num_elements;
	float target_mean;
	float target_stdev;
	float decay_constant;

} rescale_t __attribute__((__aligned__(16)));

typedef struct _rescale_gpu_t {
	/* stuff for rescaling */
	array4d_t sum;
	array4d_t sum2; // sum v**2
	array4d_t sum3; // sum v**3 - for kurtosis
	array4d_t sum4; // sum v**4 - for kurtosis
	array4d_t scale;
	array4d_t offset;
	array4d_t mean; // mean
	array4d_t std; // stdev
	array4d_t kurt; // kurtosis
	array4d_t dm0; // dm0 series
	array4d_t dm0count; // dm0 valid frequencies count
	array4d_t dm0stats; // max/min/mean/variance of Dm0 accross time
	array4d_t nsamps; // number of used samples summed
	array4d_t decay_offset;
	uint64_t interval_samps;
	uint64_t sampnum;
	int num_elements;
	float target_mean;
	float target_stdev;
	float decay_constant;
	float kurt_thresh;
	float mean_thresh;
	float std_thresh;
	float dm0_thresh;
	float cell_thresh;
	int flag_grow;
	uint64_t nf;
	uint64_t nt;
	uint64_t nbeams;

} rescale_gpu_t __attribute__((__aligned__(16)));


rescale_t* rescale_allocate(rescale_t* rescale, uint64_t nelements) ;
rescale_gpu_t* rescale_allocate_gpu(rescale_gpu_t* rescale, uint64_t nbeams, uint64_t nf, uint64_t nt, bool alloc_host);
void rescale_set_scale_offset_gpu(rescale_gpu_t* rescale, float scale, float offset);

void rescale_update_scaleoffset(rescale_t* rescale);
void rescale_update_none(rescale_t* rescale, float* inx, float*outx);
void rescale_update_float(rescale_t* rescale, float* fdata, float* sampbuf, unsigned istart);
void rescale_update_float_polsum(rescale_t* rescale, float* fdata, float* sampbuf, unsigned istart);
void rescale_update_uint8(rescale_t* rescale, float* fdata, uint8_t* sampbuf, unsigned istart);
void rescale_update_uint8_polsum(rescale_t* rescale, float* fdata, uint8_t* sampbuf, unsigned istart);
void rescale_update_int8(rescale_t* rescale, float* in, int8_t*  out);
void rescale_update_decay_float(rescale_t* rescale, float* in, float* out);
float rescale_update_decay_float_single(rescale_t* rescale, uint64_t i, float in);
void rescale_update_decay_uint8(rescale_t* rescale, float* in, uint8_t* out);
void rescale_update_scaleoffset_gpu(rescale_gpu_t& rescale);
void rescale_update_and_transpose_float_gpu(rescale_gpu_t& rescale, array4d_t& rescale_buf,
		const uint8_t* read_buf, bool invert_freq, bool subtract_dm0);

template <bool subtract_dm0> __global__ void rescale_update_and_transpose_float_kernel (
		const uint8_t* __restrict__ inarr,
		rescale_dtype* __restrict__ sumarr,
		rescale_dtype* __restrict__ sum2arr,
		rescale_dtype* __restrict__ sum3arr,
		rescale_dtype* __restrict__ sum4arr,
		rescale_dtype* __restrict__ decay_offsetarr,
		rescale_dtype* __restrict__ nsampsarr,
		const rescale_dtype* __restrict__ offsetarr,
		const rescale_dtype* __restrict__ scalearr,
		const rescale_dtype* __restrict__ dm0arr,
		const rescale_dtype* __restrict__ dm0countarr,
		const rescale_dtype* __restrict__ dm0statarr,
		rescale_dtype* __restrict__ outarr,
		float decay_constant,
		float dm0_thresh,
		float cell_thresh,
		int nt,
		bool invert_freq)
{
	int ibeam = blockIdx.x;
	int c = threadIdx.x;
	int nf = blockDim.x;
	const rescale_dtype k = decay_constant;


	// input = BTF order
	// output = BFT order
	// Rescale order: BF
	// dm0 order: BT
	// dm0sum order: B
	// nsamps order: BF

	int rsidx = c + nf*ibeam; // rescale index: BF order
	// all these reads are nice and coalesced
	rescale_dtype sum = sumarr[rsidx]; // read from global memory
	rescale_dtype sum2 = sum2arr[rsidx]; // read from global
	rescale_dtype sum3 = sum3arr[rsidx];
	rescale_dtype sum4 = sum4arr[rsidx];
	rescale_dtype decay_offset = decay_offsetarr[rsidx];  // read from global
	rescale_dtype offset = offsetarr[rsidx]; // read from global
	rescale_dtype scale = scalearr[rsidx]; // read from global
	rescale_dtype nsamps = nsampsarr[rsidx]; // read from global

	int outc;
	if (invert_freq) {
		outc = nf - 1 - c;
	} else {
		outc = c;
	}


	// Easy way of expanding the time flagging by 1. Useful for killing dropouts. ACES-209
	bool last_sample_ok = true;
	float block_dm0thresh = dm0_thresh/sqrtf((float) nt);
	int stati = 4*ibeam;
	rescale_dtype dm0min = dm0statarr[stati + 1]; // broadcast read. This is to catch dropouts
	rescale_dtype dm0stat_mean = dm0statarr[stati + 2]; // broadcast read.

	for (int t = 0; t < nt; ++t) {
		int inidx = c + nf*(t + nt*ibeam);
		int outidx = t + nt*(outc + nf*ibeam);
		// coalesced read from global
		rescale_dtype vin = (rescale_dtype)inarr[inidx]; // read from global
		rescale_dtype vout = (vin + offset) * scale;
		if (k == 0) { // If we set the timescale to zero, we just don't do any decaying
			decay_offset = 0;
		} else {
			decay_offset = (vout + decay_offset*k)/(1.0 + k);
		}
		rescale_dtype sout = vout - decay_offset;
		int dm0idx = t + nt*ibeam; // DM0 idx: BT order
		rescale_dtype dm0count = dm0countarr[dm0idx];
		rescale_dtype dm0sum = dm0arr[dm0idx] ; // sum accros dm0 - not normalised
		rescale_dtype dm0z = dm0sum*rsqrtf(dm0count);

		// the mean accross the dm0 trace is the dm0sum/dm0count (that's a mean)
		// We also subtract off the mean of the dm0 trace - the dm0state_mean has *already*
		// been normalised by rsqrtf(dm0count) - so we need todo that again to get it into a mean.
		// Bleah - I should have thorught about this harder.

		if (subtract_dm0 && scale != 0) {
			rescale_dtype dm0mean = dm0sum/dm0count - dm0stat_mean*rsqrtf(dm0count);
			sout -= dm0mean;
		}
		//int this_sample_ok = fabs(dm0) < dm0_thresh && fabs(sout) < cell_thresh && fabs(dm0sum) < block_dm0thresh;
		bool this_sample_ok = fabs(dm0z) < dm0_thresh && fabs(sout) < cell_thresh && dm0min > -3*dm0_thresh;
		//int this_sample_ok = fabs(dm0) < dm0_thresh && fabs(sout) < cell_thresh;
		if (this_sample_ok && last_sample_ok) {
			sum += vin;
			sum2 += vin*vin;
			sum3 += vin*vin*vin;
			sum4 += vin*vin*vin*vin;
			// non-coalesced write (transpose. Sorry)

			outarr[outidx] = sout;
			nsamps += 1;
		} else {
//			printf("FLAG ibeam/c/t %d/%d/%d dm0/sout/dm0min %f/%f/%f flags %d/%d/%d\n", ibeam, c, t,
//					fabs(dm0z), fabs(sout), dm0min,
//					fabs(dm0z) < dm0_thresh,
//					fabs(sout) < cell_thresh,
//					dm0min > -3*dm0_thresh);
			outarr[outidx] = 0.0;
		}

		last_sample_ok = this_sample_ok;

	}

	// write everything back to global memory -- all coalesced
	sumarr[rsidx] = sum;
	sum2arr[rsidx] = sum2;
	sum3arr[rsidx] = sum3;
	sum4arr[rsidx] = sum4;
	decay_offsetarr[rsidx] = decay_offset;
	nsampsarr[rsidx] = (float)nsamps;
}


#endif
