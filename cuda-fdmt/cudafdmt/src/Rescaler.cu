/*
 * Rescaler.cpp
 *
 *  Created on: 21 Mar 2018
 *      Author: ban115
 */

#include "Rescaler.h"

Rescaler::Rescaler(RescaleOptions& _options) : options(_options)
{
	bool alloc_host = true; // Need host memory allocated for rescale because we copy back to count flags
	int nbeams = options.nbeams;
	int nf = options.nf;
	int nt = options.nt;
	num_elements = nbeams*nf;
	sampnum = 0;
	options.interval_samps = nt;
	rescale_arraymalloc(&sum, nbeams, nf, alloc_host);
	rescale_arraymalloc(&sum2, nbeams, nf, alloc_host);
	rescale_arraymalloc(&sum3, nbeams, nf, alloc_host);
	rescale_arraymalloc(&sum4, nbeams, nf, alloc_host);
	rescale_arraymalloc(&mean, nbeams, nf, alloc_host);
	rescale_arraymalloc(&std, nbeams, nf, alloc_host);
	rescale_arraymalloc(&kurt, nbeams, nf, alloc_host);
	rescale_arraymalloc(&dm0, nbeams, nt, alloc_host);
	rescale_arraymalloc(&dm0count, nbeams, nt, alloc_host);
	rescale_arraymalloc(&dm0stats, nbeams, 4, alloc_host); // max, min, mean, var
	rescale_arraymalloc(&nsamps, nbeams, nf, alloc_host);
	rescale_arraymalloc(&scale, nbeams, nf, alloc_host);
	rescale_arraymalloc(&offset, nbeams, nf, alloc_host);
	rescale_arraymalloc(&decay_offset, nbeams, nf, alloc_host);
	array4d_set(&scale, 1.0);

}

Rescaler::~Rescaler() {

}

void Rescaler::update_scaleoffset() {
	assert(options.interval_samps > 0);
	int nf = options.nf;
	int nthreads = nf;
	assert(num_elements % nthreads == 0);
	int nblocks = num_elements / nthreads;
	rescale_update_scaleoffset_kernel<<<nblocks, nthreads>>>(
			sum.d_device,
			sum2.d_device,
			sum3.d_device,
			sum4.d_device,
			mean.d_device,
			std.d_device,
			kurt.d_device,
			offset.d_device,
			scale.d_device,
			nsamps.d_device,
			options.target_stdev,
			options.target_mean,
			options.mean_thresh,
			options.std_thresh,
			options.kurt_thresh,
			options.flag_grow);
	gpuErrchk(cudaDeviceSynchronize());
	sampnum = 0;
}

void Rescaler::set_scaleoffset(float s_scale, float s_offset) {
	array4d_set(&scale, s_scale);
	array4d_set(&offset, s_offset);
}

void Rescaler::update_and_transpose(array4d_t& rescale_buf, void* read_buf_device) {
	int nbits = options.nbits;
	switch (nbits) {
	case 1:
		do_update_and_transpose<32, uint32_t>(rescale_buf, (uint32_t*)read_buf_device);
		break;

	case 2:
		do_update_and_transpose<16, uint32_t>(rescale_buf, (uint32_t*)read_buf_device);
		break;

	case 4:
		do_update_and_transpose<8, uint32_t>(rescale_buf, (uint32_t*)read_buf_device);
		break;

	case 8:
		do_update_and_transpose<4, uint32_t>(rescale_buf, (uint32_t*)read_buf_device);
		break;

	case 16:
		do_update_and_transpose<2, uint32_t>(rescale_buf, (uint32_t*)read_buf_device);
		break;

	case 32:
		do_update_and_transpose<1, float>(rescale_buf, (float*)read_buf_device);
		break;

	default:
		printf("Invalid nbits %d for rescaler: %d\n");
		exit(1);

	}
}
