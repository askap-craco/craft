/*
 * Rescaler.cpp
 *
 *  Created on: 21 Mar 2018
 *      Author: ban115
 */

#include "Rescaler.h"
#include <fstream> // Learning C++ good boy KB!
#include <sstream>

Rescaler::Rescaler(RescaleOptions& _options, FreddaParams& _params) :
		params(_params), options(_options) {
	bool alloc_host = true; // Need host memory allocated for rescale because we copy back to count flags
	const int nb = options.nbeams_per_ant;
	const int na = options.nants;
	const int total_nbeams = na * nb;
	const int nf = options.nf;
	const int nt = options.nt;
	num_elements = total_nbeams * nf;
	num_elements_per_ant = nb * nf;
	sampnum = 0;
	options.interval_samps = nt;
	rescale_arraymalloc(&sum, na, nb, nf, alloc_host);
	rescale_arraymalloc(&sum2, na, nb, nf, alloc_host);
	rescale_arraymalloc(&sum3, na, nb, nf, alloc_host);
	rescale_arraymalloc(&sum4, na, nb, nf, alloc_host);
	rescale_arraymalloc(&mean, na, nb, nf, alloc_host);
	rescale_arraymalloc(&std, na, nb, nf, alloc_host);
	rescale_arraymalloc(&kurt, na, nb, nf, alloc_host);
	rescale_arraymalloc(&dm0, na, nb, nt, alloc_host);
	rescale_arraymalloc(&dm0count, na, nb, nt, alloc_host); // reset on update
	rescale_arraymalloc(&dm0stats, na, nb, 4, alloc_host); // max, min, mean, var
	rescale_arraymalloc(&nsamps, na, nb, nf, alloc_host);
	rescale_arraymalloc(&scale, na, nb, nf, alloc_host);
	rescale_arraymalloc(&offset, na, nb, nf, alloc_host);
	rescale_arraymalloc(&decay_offset, na, nb, nf, alloc_host);

	weights.nw = 1;
	weights.nx = na;
	weights.ny = nb;
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
	switch (options.in_order) {
	case DataOrder::BPTF:
		break; // supported for all values of NBITS

	case DataOrder::TFBP:
		assert(options.nbits == 32); // only supported for nbits=32
		break;

	default:
		printf("Invalid rescaling input order\n");
		exit(EXIT_FAILURE);
	}
	if (params.do_dump_rescaler) {
		int rescale_tsamp_mult = params.nt * params.num_rescale_blocks;
		rescale_dumpers.push_back(
				new Array4dDumper(mean, "mean", params, rescale_tsamp_mult));
		rescale_dumpers.push_back(
				new Array4dDumper(std, "std", params, rescale_tsamp_mult));
		rescale_dumpers.push_back(
				new Array4dDumper(kurt, "kurt", params, rescale_tsamp_mult));
		rescale_dumpers.push_back(
				new Array4dDumper(scale, "scale", params, rescale_tsamp_mult,
						false)); // scale is always copied
		rescale_dumpers.push_back(
				new Array4dDumper(offset, "offset", params,
						rescale_tsamp_mult));
		rescale_dumpers.push_back(
				new Array4dDumper(nsamps, "nsamps", params, rescale_tsamp_mult,
						false)); // nsamps is always copied
		rescale_dumpers.push_back(
				new Array4dDumper(decay_offset, "decay_offset", params,
						rescale_tsamp_mult));

		int block_rsamp_mult = params.nt;
		block_dumpers.push_back(
				new Array4dDumper(dm0, "dm0", params, block_rsamp_mult));
		block_dumpers.push_back(
				new Array4dDumper(dm0count, "dm0count", params,
						block_rsamp_mult));
		block_dumpers.push_back(
				new Array4dDumper(dm0stats, "dm0stats", params,
						block_rsamp_mult));


	}
	if (params.flag_file != NULL) {
		flag_frequencies_from_file(params.flag_file);
	}
}

Rescaler::~Rescaler() {
	for (auto d : rescale_dumpers) {
		delete d;
	}
	rescale_dumpers.clear();

	for (auto d : block_dumpers) {
		delete d;
	}
	block_dumpers.clear();

}

void Rescaler::reset_output(array4d_t& rescale_buf) {
	// clear output
	array4d_cuda_memset(&rescale_buf, 0);
}

void Rescaler::update_scaleoffset(RescaleOptions& options, int iant,
		cudaStream_t stream) {
	assert(options.interval_samps > 0);
	int nf = options.nf;
	int nthreads = nf;
	assert(num_elements_per_ant % nthreads == 0);
	int nblocks = num_elements_per_ant / nthreads;
	int boff = iant * options.nbeams_per_ant;
	rescale_update_scaleoffset_kernel<<<nblocks, nthreads, 0, stream>>>(
			sum.d_device, sum2.d_device, sum3.d_device, sum4.d_device,
			decay_offset.d_device, mean.d_device, std.d_device, kurt.d_device,
			offset.d_device, scale.d_device, weights.d_device, nsamps.d_device,
			options.target_stdev, options.target_mean, options.mean_thresh,
			options.std_thresh, options.kurt_thresh, options.flag_grow, boff);
	gpuErrchk(cudaPeekAtLastError());
	if (iant == 0) {
		sampnum = 0;
	}
}

void Rescaler::dump_rescale_data() {
	tdump.start();
	for (auto d : rescale_dumpers) {
		d->dump();
	}
	tdump.stop();
}

void Rescaler::dump_block_data() {
	tdump.start();
	for (auto d : block_dumpers) {
		d->dump();
	}
	tdump.stop();
}

void Rescaler::reset_ant_stats_for_first_block(int iant) {
	// reset a single antenna - dont dump or count flags
	// This is a bit hacky, but it'll do fo rnow
	assert(iant >= 0 && iant < options.nants);

	// Set nsamps array for this antenna only to zero
	size_t size =  nsamps.ny * nsamps.nz;
	size_t idx = array4d_idx(&nsamps, 0, iant, 0, 0);
	assert(nsamps.d_device != NULL);
	gpuErrchk(cudaMemset(&nsamps.d_device[idx], 0, size * sizeof(fdmt_dtype)));

	// Copy block mean to rescaler mean so it starts near the right place when we decode this block
	gpuErrchk(cudaMemcpy(&decay_offset.d_device[idx], &mean.d_device[idx], size*sizeof(fdmt_dtype), cudaMemcpyDeviceToDevice));
}

void Rescaler::finish_all_ants() {
	// needs to be called once all blocks have finished. e.g. after a cudaDeviceSynchronize()

	if (params.do_dump_rescaler) {
		dump_block_data(); // Run this every block
	}
}

void Rescaler::update_rescale_statistics() {
	for (int iant = 0; iant < params.source->nants(); ++iant) {
		update_scaleoffset(options, iant);
	}

	// Count how many  channels have been flagged for this whole block
	// by looking at how many channels have scale==0
	array4d_copy_to_host(&scale);
	for (int i = 0; i < params.nf * params.nbeams_in_total; ++i) {
		if (scale.d[i] == 0) {
			// that channel will stay flagged for num_rescale_blocks
			num_flagged_beam_chans += params.num_rescale_blocks;
		}
	}

	// Count how many times were flagged
	//assert(params.num_rescale_blocks >= 0);
	array4d_copy_to_host(&nsamps); // must do this before updating scaleoffset, which resets nsamps to zero

	// Calculate number of DM0 times that were flagged
	for (int i = 0; i < params.nf * params.nbeams_in_total; ++i) {
		int ns = (int) nsamps.d[i]; // nsamps is the number of unflagged samples from this block
		int nflagged = params.nt*params.num_rescale_blocks - ns;
		// rescale.sampnum is the total number of samples that has gone into the rescaler since resetting
		assert(nflagged >= 0);
		num_flagged_times += nflagged;
	}


	if (params.do_dump_rescaler) {
		dump_rescale_data();
	}

	// reset nsamps
	array4d_cuda_memset(&nsamps, 0);

}

void Rescaler::set_scaleoffset(float s_scale, float s_offset) {
	array4d_set(&scale, s_scale);
	array4d_set(&offset, s_offset);
}

int Rescaler::flag_frequencies_from_file(const char* filename) {
	std::ifstream infile(filename);
	std::string line;
	float freq;
	int num_flagged = 0;
	if (infile.is_open()) {
		while (std::getline(infile, line)) {
			if (line[0] == '#') {
				continue;
			}
			std::istringstream iss(line);
			if (iss >> freq) {
				bool was_flagged = flag_frequency(freq);
				if (was_flagged) {
					num_flagged += 1;
				}
			}
		}
	} else {
		perror("Could not open flag file file");
		exit(EXIT_FAILURE);
	}

	printf("Flagged %d channels from flagfile %s\n", num_flagged, filename);

	return num_flagged;
}

bool Rescaler::flag_frequency(float freq) {
	const float fch1 = params.source->fch1();
	const float foff = params.source->foff();
	int channel = int(roundf((freq - fch1) / foff));
	bool flagged = false;
	if (channel >= 0 && channel < params.nf) {
		printf("Flagging channel %d at frequency %f\n", channel, freq);
		flag_channel(channel);
		flagged = true;
	}

	return flagged;
}

void Rescaler::flag_channel(int channel) {
	assert(channel >= 0 && channel < options.nf);
	for (int iant = 0; iant < options.nants; iant++) {
		for (int ibeam = 0; ibeam < options.nbeams_per_ant; ibeam++) {
			int idx = array4d_idx(&weights, 0, iant, ibeam, channel);
			weights.d[idx] = 0.0f;
		}
	}
	array4d_copy_to_device(&weights);
}

void Rescaler::flag_beam(int beam) {
	assert(beam >=0 && beam < options.nbeams_per_ant);
	for (int iant = 0; iant < options.nants; iant++) {
		for(int ichan = 0; ichan < options.nf; ichan++) {
			int idx = array4d_idx(&weights, 0, iant, beam, ichan);
			weights.d[idx] = 0.0f;
		}
	}
	array4d_copy_to_device(&weights);
}

void Rescaler::process_ant_block(array4d_t& rescale_buf, void* read_buf_device,
		RescaleOptions &options, int iant, cudaStream_t stream) {
	int nbits = options.nbits;
	switch (nbits) {
	case 1:
		do_update_and_transpose<32, uint32_t>(rescale_buf,
				(uint32_t*) read_buf_device, options, iant, stream);
		break;

	case 2:
		do_update_and_transpose<16, uint32_t>(rescale_buf,
				(uint32_t*) read_buf_device, options, iant, stream);
		break;

	case 4:
		do_update_and_transpose<8, uint32_t>(rescale_buf,
				(uint32_t*) read_buf_device, options, iant, stream);
		break;

	case 8:
		do_update_and_transpose<4, uint32_t>(rescale_buf,
				(uint32_t*) read_buf_device, options, iant, stream);
		break;

	case 16:
		do_update_and_transpose<2, uint32_t>(rescale_buf,
				(uint32_t*) read_buf_device, options, iant, stream);
		break;

	case 32:
		do_update_and_transpose<1, float>(rescale_buf, (float*) read_buf_device,
				options, iant, stream);
		break;

	default:
		printf("Invalid nbits %d for rescaler: %d\n");
		exit(1);

	}
}
