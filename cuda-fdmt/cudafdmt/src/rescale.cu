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
		float meansq = rescale->sum2[i]/nsamp;
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
		rescale->sum2[i] = 0.0;
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
		rescale->sum2[i] += vin2;
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
		rescale->sum2[j] += vin*vin;

		float uin = inx[j+1];
		rescale->sum[j+1] += uin;
		rescale->sum2[j+1] += uin*uin;

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
		rescale->sum2[i] += vin*vin;
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
		rescale->sum2[j] += vin*vin;

		float uin=in[j+1];
		rescale->sum[j+1] += uin;
		rescale->sum2[j+1] += uin*uin;

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
		rescale->sum2[i] += vin*vin;
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
	rescale->sum2[i] += vin*vin;
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
		rescale->sum2[i] += vin*vin;
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
		rescale->sum2[i] += vin*vin;
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
	rescale->sum2 = (float*) rescale_malloc(sz);
	rescale->scale = (float*) rescale_malloc(sz);
	rescale->offset = (float*) rescale_malloc(sz);
	rescale->decay_offset = (float*) rescale_malloc(sz);
	rescale->sampnum = 0;
	rescale->num_elements = nelements;

	for(uint64_t i = 0; i < nelements; ++i) {
		rescale->sum[i] = 0;
		rescale->sum2[i] = 0;
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

void rescale_arraymalloc(array4d_t* arr, uint64_t nbeams, uint64_t nf, bool alloc_host)
{
	arr->nw = 1;
	arr->nx = 1;
	arr->ny = nbeams;
	arr->nz = nf;
	array4d_malloc(arr, alloc_host, true);
	array4d_zero(arr);
}

rescale_gpu_t* rescale_allocate_gpu(rescale_gpu_t* rescale, uint64_t nbeams, uint64_t nf, uint64_t nt, bool alloc_host)
{
	uint64_t nelements = nbeams*nf;
	size_t sz = nelements*sizeof(rescale_dtype);
	rescale_arraymalloc(&rescale->sum, nbeams, nf, alloc_host);
	rescale_arraymalloc(&rescale->sum2, nbeams, nf, alloc_host);
	rescale_arraymalloc(&rescale->sum3, nbeams, nf, alloc_host);
	rescale_arraymalloc(&rescale->sum4, nbeams, nf, alloc_host);
	rescale_arraymalloc(&rescale->mean, nbeams, nf, alloc_host);
	rescale_arraymalloc(&rescale->std, nbeams, nf, alloc_host);
	rescale_arraymalloc(&rescale->kurt, nbeams, nf, alloc_host);
	rescale_arraymalloc(&rescale->dm0, nbeams, nt, alloc_host);
	rescale_arraymalloc(&rescale->dm0count, nbeams, nt, alloc_host);
	rescale_arraymalloc(&rescale->dm0stats, nbeams, 4, alloc_host); // max, min, mean, var
	rescale_arraymalloc(&rescale->nsamps, nbeams, nf, alloc_host);
	rescale_arraymalloc(&rescale->scale, nbeams, nf, alloc_host);
	rescale_arraymalloc(&rescale->offset, nbeams, nf, alloc_host);
	rescale_arraymalloc(&rescale->decay_offset, nbeams, nf, alloc_host);
	array4d_set(&rescale->scale, 1.0);

	rescale->sampnum = 0;
	rescale->num_elements = nelements;
	rescale->nf = nf;
	rescale->nt = nt;
	rescale->nbeams = nbeams;

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


__global__ void rescale_calc_dm0_kernel (
		const uint8_t* __restrict__ inarr,
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
	for(int t = threadIdx.x; t < nt; t += blockDim.x) {
		rescale_dtype dm0sum = 0.0;
		int nsamp = 0;
		for (int c = 0; c < nf; ++c) {
			int rsidx = c + nf*ibeam; // rescale index BF order
			// all these reads are nice and coalesced
			rescale_dtype offset = offsetarr[rsidx]; // read from global
			rescale_dtype scale = scalearr[rsidx]; // read from global
			int inidx = c + nf*(t + nt*ibeam); // input index : BTF order

			// coalesced read from global
			rescale_dtype vin = (rescale_dtype)inarr[inidx]; // read from global
			rescale_dtype vout = (vin + offset) * scale;
			if (fabs(vout) < cell_thresh && scale != 0.0f) {
				dm0sum += vout;
				++nsamp;
			}
		}

		int dm0idx = t + nt*ibeam;
		rescale_dtype correction = rsqrtf((float) nsamp);
		//dm0arr[dm0idx] = dm0sum * correction;
		dm0arr[dm0idx] = dm0sum;
		dm0count[dm0idx] = (float)nsamp;
	}
}

__global__ void rescale_calc_dm0stats_kernel (
		const rescale_dtype* __restrict__ dm0arr,
		const rescale_dtype* __restrict__ dm0countarr,
		rescale_dtype* __restrict__ dm0statarr,
		int nt)
{
	// dm0 order: BT
	// dm0stats order: BX
	// X is max,min,mean,var

	int ibeam = threadIdx.x;
	rescale_dtype dm0sum = 0.0;
	rescale_dtype dm0sum2 = 0.0;
	rescale_dtype nsampinit = dm0countarr[ibeam];
	rescale_dtype vinit =  dm0arr[ibeam] * rsqrtf(nsampinit); // normalise to sqrt number of additions
	rescale_dtype dm0min = vinit;
	rescale_dtype dm0max = vinit;


	for (int t = 0; t < nt; ++t) {
		int dmidx = t + nt*ibeam;
		rescale_dtype nsamp = dm0countarr[dmidx];
		rescale_dtype v = dm0arr[dmidx] * rsqrtf(nsamp); // normalise to sqrt number of additions

		dm0sum += v;
		dm0sum2 += v*v;
		if (v < dm0min) {
			dm0min = v;
		}
		if (v > dm0max) {
			dm0max = v;
		}

	}
	//dm0sumarr[ibeam] = dm0sum/((float) nt);
	rescale_dtype nsamp = (float) nt;
	rescale_dtype dm0mean = dm0sum/nsamp;
	rescale_dtype mean2 = dm0sum2/nsamp;
	rescale_dtype dm0var = mean2 - dm0mean*dm0mean;

	dm0statarr[ibeam + 0] = dm0max;
	dm0statarr[ibeam + 1] = dm0min;
	dm0statarr[ibeam + 2] = dm0mean;
	dm0statarr[ibeam + 3] = dm0var;

	//printf("DM stats ibeam=%d max/min/mean/var %f/%f/%f/%f\n", ibeam, dm0max, dm0min, dm0mean, dm0var);
}


__global__ void rescale_update_and_transpose_float_kernel (
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
		bool invert_freq,
		bool subtract_dm0)
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
	rescale_dtype dm0min = dm0statarr[ibeam + 1]; // broadcast read. This is to catch dropouts

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
		rescale_dtype dm0mean = dm0sum/dm0count;
		//int this_sample_ok = fabs(dm0) < dm0_thresh && fabs(sout) < cell_thresh && fabs(dm0sum) < block_dm0thresh;
		bool this_sample_ok = fabs(dm0z) < dm0_thresh && fabs(sout) < cell_thresh && dm0min > -3*dm0_thresh;
		//int this_sample_ok = fabs(dm0) < dm0_thresh && fabs(sout) < cell_thresh;
		if (this_sample_ok && last_sample_ok) {
			sum += vin;
			sum2 += vin*vin;
			sum3 += vin*vin*vin;
			sum4 += vin*vin*vin*vin;
			// non-coalesced write (transpose. Sorry)
			if (subtract_dm0) {
				outarr[outidx] = sout - dm0mean;
			} else {
				outarr[outidx] = sout;
			}
			nsamps += 1;
		} else {
			printf("NOK ibeam/c/t %d/%d/%d dm0/sout/dm0min %f/%f/%f flags %d/%d/%d\n", ibeam, c, t,
					fabs(dm0z), fabs(sout), dm0min,
					fabs(dm0z) < dm0_thresh,
					fabs(sout) < cell_thresh,
					dm0min > -3*dm0_thresh);
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

void rescale_update_and_transpose_float_gpu(rescale_gpu_t& rescale, array4d_t& rescale_buf, const uint8_t* read_buf, bool invert_freq, bool subtract_dm0)
{
	int nbeams = rescale_buf.nw;
	int nf = rescale_buf.nx;
	int nt = rescale_buf.nz;

	// Calculate dm0 for flagging
	rescale_calc_dm0_kernel<<<nbeams, 256>>>(
			read_buf,
			rescale.offset.d_device,
			rescale.scale.d_device,
			rescale.dm0.d_device,
			rescale.dm0count.d_device,
			nf, nt,
			rescale.cell_thresh);

	gpuErrchk(cudaDeviceSynchronize());

	// Take the mean all the dm0 times into one big number per beam - this is the how we flag
	// short dropouts see ACES-209
	// probably could do this in rescale_calc_dm0_kernel after yu've done it
	// But i Haven't got htere yet.
	rescale_calc_dm0stats_kernel<<<1, nbeams>>>(
			rescale.dm0.d_device,
			rescale.dm0count.d_device,
			rescale.dm0stats.d_device,
			nt);

	gpuErrchk(cudaDeviceSynchronize());

	rescale_update_and_transpose_float_kernel<<<nbeams, nf>>>(
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
		int flag_grow)
{
	int c = threadIdx.x;
	int nf = blockDim.x;
	int ibeam = blockIdx.x;
	int i = c + nf*ibeam;
	rescale_dtype nsamp = nsamparr[i];
	rescale_dtype mean = sum[i]/nsamp;
	rescale_dtype mean2 = sum2[i]/nsamp;
	rescale_dtype mean3 = sum3[i]/nsamp;
	rescale_dtype mean4 = sum4[i]/nsamp;
	rescale_dtype variance = mean2 - mean*mean;

	// Excess Kurtosis is k = E([X-mu]**4)/(Var[X]**2) - 3
	// numerator = E[X**4] - 4E[X][E[X**3] + 6 E[X**2]E[X]**2 - 3E[X]**4
	rescale_dtype kurt = (mean4 - 4*mean*mean3 + 6*mean2*mean*mean - 3*mean*mean*mean*mean)/(variance*variance) -3 ;
	// save flag inputs
	meanarr[i] = mean;
	stdarr[i] = sqrtf(variance);
	kurtarr[i] = kurt;

	__syncthreads();
	rescale_dtype scale = 0.0, offset = 0.0;
	int icstart = max(0, c - flag_grow) + nf*ibeam;
	int icend = min(nf, c + flag_grow) + nf*ibeam;
	int flag = 0;
	rescale_dtype expected_mean = 128.;
	rescale_dtype expected_std = 18.;
	for (int ic = icstart; ic < icend; ++ic) {
		rescale_dtype meanoff = fabs(meanarr[ic] - expected_mean);
		rescale_dtype stdoff = fabs(stdarr[ic] - expected_std);
		rescale_dtype kurtoff = fabs(kurtarr[ic]);

		if (nsamp > 0 && (meanoff > mean_thresh ||
				stdoff > std_thresh ||
				kurtoff > kurt_thresh)) {
			flag = 1;
			break;
		}
	}

	if (flag) {
		scale = 0.0;
		offset = 0.0;
	} else {
		if (nsamp == 0) { // Don't update the scale and offset if everything has been flagged
			offset = offsetarr[i];
			scale = scalearr[i];
		} else if (variance == 0.0) {
			scale = 1.0;
			offset = -mean + target_mean/scale;
		} else {
			scale = target_stdev / sqrt(variance);
			offset = -mean + target_mean/scale;
		}
	}

	offsetarr[i] = offset;
	scalearr[i] = scale;

	// reset values to zero
	sum[i] = 0.0;
	sum2[i] = 0.0;
	sum3[i] = 0.0;
	sum4[i] = 0.0;
	nsamparr[i] = 0;
}
void rescale_update_scaleoffset_gpu(rescale_gpu_t& rescale)
{
	assert(rescale.interval_samps > 0);
	int nthreads = rescale.nf;
	assert(rescale.num_elements % nthreads == 0);
	int nblocks = rescale.num_elements / nthreads;
	rescale_update_scaleoffset_kernel<<<nblocks, nthreads>>>(
			rescale.sum.d_device,
			rescale.sum2.d_device,
			rescale.sum3.d_device,
			rescale.sum4.d_device,
			rescale.mean.d_device,
			rescale.std.d_device,
			rescale.kurt.d_device,
			rescale.offset.d_device,
			rescale.scale.d_device,
			rescale.nsamps.d_device,
			rescale.target_stdev,
			rescale.target_mean,
			rescale.mean_thresh,
			rescale.std_thresh,
			rescale.kurt_thresh,
			rescale.flag_grow);
	gpuErrchk(cudaDeviceSynchronize());
	rescale.sampnum = 0;
}

