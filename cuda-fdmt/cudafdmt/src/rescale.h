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

void rescale_arraymalloc(array4d_t* arr, uint64_t nbeams, uint64_t nf, bool alloc_host);

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
__global__ void rescale_calc_dm0stats_kernel (
		const rescale_dtype* __restrict__ dm0arr,
		const rescale_dtype* __restrict__ dm0countarr,
		rescale_dtype* __restrict__ dm0statarr,
		int nt);

__global__ void rescale_update_scaleoffset_kernel (
		rescale_dtype* __restrict__ sum,
		rescale_dtype* __restrict__ sum2,
		rescale_dtype* __restrict__ sum3,
		rescale_dtype* __restrict__ sum4,
		rescale_dtype* __restrict__ meanarr,
		rescale_dtype* __restrict__ stdarr,
		rescale_dtype* __restrict__ kurtarr,
		rescale_dtype* __restrict__ offsetarr,
		rescale_dtype* __restrict__ scalearr,
		rescale_dtype* nsamparr,
		rescale_dtype target_stdev,
		rescale_dtype target_mean,
		rescale_dtype mean_thresh,
		rescale_dtype std_thresh,
		rescale_dtype kurt_thresh,
		int flag_grow);

template <int nsamps_per_word, typename wordT> __device__ __host__ inline rescale_dtype extract_sample(const wordT word, const int samp)
{
	const int nbits_per_samp = sizeof(wordT)*8/nsamps_per_word;
	const int shift = nbits_per_samp * samp;

	// Shift desired sample down to bottom of shiftsamp. Most significant bits are sample=0
	wordT shiftsamp = word >> shift;

	// this mask has 1s in the bottom nbits_per_samp bits
	wordT mask = (1 << nbits_per_samp) - 1;

	wordT masksamp = shiftsamp & mask;

	rescale_dtype sample = (rescale_dtype) masksamp;
	/*
	printf("nbits_per_samp %d shift %d shiftsamp %d nsamp_per_word %d word %d mask %d masksamp %d sample %f\n",
			nbits_per_samp, shift, shiftsamp, nsamps_per_word, word, mask, masksamp, sample);
			*/


	return sample;

}

template <> __device__ __host__ inline rescale_dtype extract_sample<1, float>(const float word, int sampno)
{
	return (rescale_dtype) word;
}

template <> __device__ __host__ inline rescale_dtype extract_sample<1, double>(const double word, int sampno)
{
	return (rescale_dtype) word;
}

template <int nsamps_per_word, typename wordT> __global__ void rescale_update_and_transpose_float_kernel (
		const wordT* __restrict__ inarr,
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
		bool invert_freq,
		bool subtract_dm0,
		bool polsum)
{
	int ibeam = blockIdx.x;
	int s = threadIdx.x; // sample index within a word
	int w = threadIdx.y; // word index
	int nwords = blockDim.y;
	const int nf = nwords * nsamps_per_word;
	const int c = w*nsamps_per_word + s; // channel number
	const rescale_dtype k = decay_constant;

	int outbeam;
	if (polsum) {
		// TODO: Maybe always polsum and do npol=1 and beams*2 or whatever
		outbeam = ibeam/2; // ASSUMES BP  order and npol = 2
	} else {
		outbeam = ibeam;
	}


	// on input F axis is broken further into words and samples. all X threads load the same word
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
		int wordidx = w + nwords*(t + nt*ibeam);

		// coalesced read from global for all x threads.
		wordT word = inarr[wordidx];
		rescale_dtype vin = extract_sample<nsamps_per_word, wordT>(word, s);

		rescale_dtype vout = (vin + offset) * scale;
		rescale_dtype sout = vout;
		if (k == 0) { // If we set the timescale to zero, we just don't do any decaying
			decay_offset = 0;
			sout = vout;
		} else {
			decay_offset = (vout + decay_offset*k)/(1.0 + k);
			sout = vout - decay_offset;
		}
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
		int outidx = t + nt*(outc + nf*outbeam);

		if (this_sample_ok && last_sample_ok) {
			sum += vin;
			sum2 += vin*vin;
			sum3 += vin*vin*vin;
			sum4 += vin*vin*vin*vin;
			// non-coalesced write (transpose. Sorry)

			//outarr[outidx] += sout;

			// doing an atomic add for polarisation summing and antenna summing - probably adds some overhead but we'll see.t
			atomicAdd(outarr + outidx, sout);
			nsamps += 1;
		} else {
//			printf("FLAG ibeam/c/t %d/%d/%d dm0/sout/dm0min %f/%f/%f flags %d/%d/%d\n", ibeam, c, t,
//					fabs(dm0z), fabs(sout), dm0min,
//					fabs(dm0z) < dm0_thresh,
//					fabs(sout) < cell_thresh,
//					dm0min > -3*dm0_thresh);
			//outarr[outidx] = 0.0;
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

template <int nsamps_per_word, typename wordT> __global__ void rescale_calc_dm0_kernel (
		const wordT* __restrict__ inarr,
		const rescale_dtype* __restrict__ offsetarr,
		const rescale_dtype* __restrict__ scalearr,
		rescale_dtype* __restrict__ dm0arr,
		rescale_dtype* __restrict__ dm0count,
		int nf,
		int nt,
		rescale_dtype cell_thresh)
{
	// input = BTF order
	// dm0 order: BT
	// Rescale: BF order

	int ibeam = blockIdx.x;

	// to take advantage of having a word in a register, we have to have two loops to loop
	// over channels - one loop over the words, and the next over the samples in the word
	int nwords = nf/nsamps_per_word;

	for(int t = threadIdx.x; t < nt; t += blockDim.x) {
		rescale_dtype dm0sum = 0.0;
		int nsamp = 0;
		for (int w = 0; w < nwords; w++) {
			int inidx = w + nwords*(t + nt*ibeam); // input index : BTF order
			// coalesced read from global
			wordT word = inarr[inidx];

			for (int s = 0; s < nsamps_per_word; ++s) {
				int c = w*nwords + s; // channel number
				int rsidx = c + nf*ibeam; // rescale index BF order
				// all these reads are nice and coalesced
				rescale_dtype offset = offsetarr[rsidx]; // read from global
				rescale_dtype scale = scalearr[rsidx]; // read from global

				// extract channel out of word
				rescale_dtype vin = extract_sample<nsamps_per_word, wordT>(word, s); // read from global
				rescale_dtype vout = (vin + offset) * scale;
				if (fabs(vout) < cell_thresh && scale != 0.0f) {
					dm0sum += vout;
					++nsamp;
				}

			}
		}

		int dm0idx = t + nt*ibeam;
		rescale_dtype correction = rsqrtf((float) nsamp);
		//dm0arr[dm0idx] = dm0sum * correction;
		dm0arr[dm0idx] = dm0sum;
		dm0count[dm0idx] = (float)nsamp;
	}
}



template <int nsamps_per_word, typename wordT> void
	rescale_update_and_transpose_float_gpu(rescale_gpu_t& rescale,
										array4d_t& rescale_buf,
										const wordT* read_buf,
										bool invert_freq,
										bool subtract_dm0)
{
	int nbeams = rescale_buf.nw;
	int nf = rescale_buf.nx;
	int nt = rescale_buf.nz;
	int nwords = nf / nsamps_per_word;
	assert(nf % nsamps_per_word == 0);

	// clear output
	array4d_cuda_memset(&rescale_buf, 0);

	rescale_calc_dm0_kernel< nsamps_per_word, wordT > <<<nbeams, 256>>>(
			read_buf,
			rescale.offset.d_device,
			rescale.scale.d_device,
			rescale.dm0.d_device,
			rescale.dm0count.d_device,
			nf, nt,
			rescale.cell_thresh);

	// Take the mean all the dm0 times into one big number per beam - this is the how we flag
	// short dropouts see ACES-209
	// probably could do this in rescale_calc_dm0_kernel after yu've done it
	// But i Haven't got htere yet.
	rescale_calc_dm0stats_kernel<<<1, nbeams>>>(
			rescale.dm0.d_device,
			rescale.dm0count.d_device,
			rescale.dm0stats.d_device,
			nt);

	dim3 blockdim(nsamps_per_word, nwords);

	rescale_update_and_transpose_float_kernel< nsamps_per_word, wordT ><<<nbeams, blockdim>>>(
			read_buf,
			rescale.sum.d_device,
			rescale.sum2.d_device,
			rescale.sum3.d_device,
			rescale.sum4.d_device,
			rescale.decay_offset.d_device,
			rescale.nsamps.d_device,
			rescale.offset.d_device,
			rescale.scale.d_device,
			rescale.dm0.d_device,
			rescale.dm0count.d_device,
			rescale.dm0stats.d_device,
			rescale_buf.d_device,
			rescale.decay_constant,
			rescale.dm0_thresh,
			rescale.cell_thresh*rescale.target_stdev,
			nt,
			invert_freq,
			subtract_dm0);


	rescale.sampnum += nt;
	gpuErrchk(cudaDeviceSynchronize());

}

#endif
