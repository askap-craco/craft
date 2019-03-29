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
	//printf("Allocating array ");
	//array4d_print_shape(arr);
	//printf("for rescaling\n");
	array4d_malloc(arr, alloc_host, true);
	array4d_zero(arr);
}

rescale_gpu_t* rescale_allocate_gpu(rescale_gpu_t* rescale, uint64_t nbeams, uint64_t nf, uint64_t nt, bool alloc_host)
{
	uint64_t nelements = nbeams*nf;
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

__global__ void rescale_calc_dm0stats_kernel (
		const rescale_dtype* __restrict__ dm0arr,
		const rescale_dtype* __restrict__ dm0countarr,
		rescale_dtype* __restrict__ dm0statarr,
		int nt,
		int boff)
{
	// dm0 order: BT
	// dm0stats order: BX
	// X is max,min,mean,var

	int ibeam = threadIdx.x;
	int rsbeam = ibeam + boff;
	rescale_dtype dm0sum = 0.0;
	rescale_dtype dm0sum2 = 0.0;
	rescale_dtype nsampinit = dm0countarr[rsbeam*nt];
	rescale_dtype vinit =  dm0arr[rsbeam*nt] * rsqrtf(nsampinit); // normalise to sqrt number of additions
	rescale_dtype dm0min = vinit;
	rescale_dtype dm0max = vinit;


	for (int t = 0; t < nt; ++t) {
		int dmidx = t + nt*rsbeam;
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

	int stati = rsbeam*4;
	dm0statarr[stati + 0] = dm0max;
	dm0statarr[stati + 1] = dm0min;
	dm0statarr[stati + 2] = dm0mean;
	dm0statarr[stati + 3] = dm0var;

	//printf("DM stats ibeam=%d max/min/mean/var %f/%f/%f/%f\n", ibeam, dm0max, dm0min, dm0mean, dm0var);
}




__global__ void rescale_update_scaleoffset_kernel (
		rescale_dtype* __restrict__ m1,
		rescale_dtype* __restrict__ m2,
		rescale_dtype* __restrict__ m3,
		rescale_dtype* __restrict__ m4,
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
		int flag_grow,
		int boff)
{
	int c = threadIdx.x;
	int nf = blockDim.x;
	int ibeam = blockIdx.x;
	int rsbeam = ibeam + boff;
	int i = c + nf*rsbeam;
	rescale_dtype nsamp = nsamparr[i];
	rescale_dtype mean = m1[i];
//	rescale_dtype mean2 = sum2[i]/nsamp;
//	rescale_dtype mean3 = sum3[i]/nsamp;
//	rescale_dtype mean4 = sum4[i]/nsamp;
	rescale_dtype variance = m2[i]/(nsamp - 1.0f);
	rescale_dtype std = sqrtf(variance);


	// Excess Kurtosis is k = E([X-mu]**4)/(Var[X]**2) - 3
	// numerator = E[X**4] - 4E[X][E[X**3] + 6 E[X**2]E[X]**2 - 3E[X]**4
	//rescale_dtype kurt = (mean4 - 4*mean*mean3 + 6*mean2*mean*mean - 3*mean*mean*mean*mean)/(variance*variance) -3
	// from here: https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance

    rescale_dtype kurt = (nsamp*m4[i]) / (m2[i]*m2[i]) - 3.0f;
	if (!isfinite(kurt)) {
		kurt = 0;
	}

	rescale_dtype skew = sqrtf(nsamp)*m3[i]/powf(m2[i], 1.5f);

	// save flag inputs

	rescale_dtype prev_mean = meanarr[i];
	rescale_dtype prev_std = stdarr[i];
	meanarr[i] = mean;
	stdarr[i] = std;
	kurtarr[i] = kurt;

	rescale_dtype scale = 0.0, offset = 0.0;
	//int icstart = max(0, c - flag_grow) + nf*rsbeam;
	//int icend = min(nf, c + flag_grow) + nf*rsbeam;
	rescale_dtype meanoff, stdoff, kurtoff;
	bool flag = false;
	if (prev_mean == 0) {
		meanoff = 0;
	} else {
		meanoff = fabs(mean - prev_mean)/prev_mean;
	}
	if (prev_std == 0) {
		stdoff = 0;
	} else {
		stdoff = fabs(std - prev_std)/prev_std;
	}
	kurtoff = fabs(kurt);

	// some of these divisions can be by 0, and the thresholds can be inf - this handles all that.

	 bool thresh_ok = (meanoff <= mean_thresh) &&
			(stdoff <= std_thresh) &&
			(kurtoff <= kurt_thresh);

	if (! thresh_ok) {
//			printf("Rescale ibeam %d ic=%d meanoff=%f prev_mean=%f meanarr=%f mean_thresh=%f stdoff=%f prev_std=%f stdarr=%f kurtoff=%f thresh OK? %d\n",
//							ibeam, ic, meanoff, prev_mean, meanarr[ic], mean_thresh, stdoff, prev_std, stdarr[ic], kurtoff, thresh_ok);

		flag = true;

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
			scale = target_stdev / sqrtf(variance);
			offset = -mean + target_mean/scale;
		}
	}

	offsetarr[i] = offset;
	scalearr[i] = scale;

	// reset values to zero
	m1[i] = 0.0;
	m2[i] = 0.0;
	m3[i] = 0.0;
	m4[i] = 0.0;
	nsamparr[i] = 0;
}
void rescale_update_scaleoffset_gpu(rescale_gpu_t& rescale, int iant)
{
	assert(rescale.interval_samps > 0);
	int nthreads = rescale.nf;
	assert(rescale.num_elements % nthreads == 0);
	int nblocks = rescale.num_elements / nthreads;
	int boff = iant*rescale.nbeams;
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
			rescale.flag_grow,
			boff);
	// zero decay offsets after updating block offsets - otherwise you get big steps. CRAFT-206
	array4d_zero(&rescale.decay_offset);
	gpuErrchk(cudaDeviceSynchronize());
	rescale.sampnum = 0;
}

