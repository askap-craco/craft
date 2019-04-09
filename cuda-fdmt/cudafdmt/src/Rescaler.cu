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
	int nbeams = options.nbeams_per_ant * options.nants;
	int nf = options.nf;
	int nt = options.nt;
	num_elements = nbeams*nf;
	num_elements_per_ant = options.nbeams_per_ant * nf;
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

	weights.nw = 1;
	weights.nx = options.nants;
	weights.ny = nbeams;
	weights.nz = nf;
	array4d_malloc(&weights, true, true);
	array4d_set(&weights, 1.0);

	array4d_set(&scale, 1.0);

	// set up the options to use if we don't want to do flagging
	noflag_options = options; // make a copy
	// For not flagging, we set the thresholds to max value
	noflag_options.mean_thresh = INFINITY;
	noflag_options.std_thresh = INFINITY;
	noflag_options.kurt_thresh = INFINITY;
	noflag_options.dm0_thresh = INFINITY;
	noflag_options.cell_thresh = INFINITY;
	noflag_options.decay_constant = 0;
	noflag_options.flag_grow = 1;
	switch(options.in_order) {
	case DataOrder::BPTF:
		break; // supported for all values of NBITS

	case DataOrder::TFBP:
		assert(options.nbits == 32); // only supported for nbits=32
		break;

	default:
		printf("Invalid rescaling input order\n");
		exit(EXIT_FAILURE);
	}
}

Rescaler::~Rescaler() {

}

void Rescaler::reset(array4d_t& rescale_buf)
{
	// clear output
	array4d_cuda_memset(&rescale_buf, 0);
}

void Rescaler::update_scaleoffset(RescaleOptions& options, int iant, cudaStream_t stream) {
	assert(options.interval_samps > 0);
	int nf = options.nf;
	int nthreads = nf;
	assert(num_elements_per_ant % nthreads == 0);
	int nblocks = num_elements_per_ant / nthreads;
	int boff = iant*options.nbeams_per_ant;
	rescale_update_scaleoffset_kernel<<<nblocks, nthreads, 0, stream>>>(
			sum.d_device,
			sum2.d_device,
			sum3.d_device,
			sum4.d_device,
			decay_offset.d_device,
			mean.d_device,
			std.d_device,
			kurt.d_device,
			offset.d_device,
			scale.d_device,
			weights.d_device,
			nsamps.d_device,
			options.target_stdev,
			options.target_mean,
			options.mean_thresh,
			options.std_thresh,
			options.kurt_thresh,
			options.flag_grow,
			boff);
	if (iant == 0) {
		sampnum = 0;
	}
}

void Rescaler::set_scaleoffset(float s_scale, float s_offset) {
	array4d_set(&scale, s_scale);
	array4d_set(&offset, s_offset);
}

void Rescaler::flag_channel(int channel) {
	assert(channel > 0 && channel < options.nf);
	for (int iant = 0 ;iant < options.nants; iant++) {
		for (int ibeam = 0; ibeam < options.nbeams_per_ant; ibeam++) {
			int idx = array4d_idx(&weights, 0, iant, ibeam, channel);
			weights.d[idx] = 0.0f;
		}
	}
	array4d_copy_to_device(&weights);
}

void Rescaler::update_and_transpose(array4d_t& rescale_buf, void* read_buf_device, RescaleOptions &options, int iant, cudaStream_t stream) {
	int nbits = options.nbits;
	switch (nbits) {
	case 1:
		do_update_and_transpose<32, uint32_t>(rescale_buf, (uint32_t*)read_buf_device, options, iant, stream);
		break;

	case 2:
		do_update_and_transpose<16, uint32_t>(rescale_buf, (uint32_t*)read_buf_device, options, iant, stream);
		break;

	case 4:
		do_update_and_transpose<8, uint32_t>(rescale_buf, (uint32_t*)read_buf_device, options, iant, stream);
		break;

	case 8:
		do_update_and_transpose<4, uint32_t>(rescale_buf, (uint32_t*)read_buf_device, options, iant, stream);
		break;

	case 16:
		do_update_and_transpose<2, uint32_t>(rescale_buf, (uint32_t*)read_buf_device, options, iant, stream);
		break;

	case 32:
		do_update_and_transpose<1, float>(rescale_buf, (float*)read_buf_device, options, iant, stream);
		break;

	default:
		printf("Invalid nbits %d for rescaler: %d\n");
		exit(1);

	}
}
