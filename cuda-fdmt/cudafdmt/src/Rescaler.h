/*
 * Rescaler.h
 *
 *  Created on: 21 Mar 2018
 *      Author: ban115
 */

#ifndef RESCALER_H_
#define RESCALER_H_
#include <vector>

#include "rescale.h"
#include "array.h"
#include "DataOrder.h"
#include "FreddaParams.h"
#include "Array4dDumper.h"
#include "cuda_utils.h"
#include "CudaTimer.h"

using std::vector;

struct RescaleOptions {
public:
	float target_mean;
	float target_stdev;
	float decay_constant;
	float mean_thresh;
	float std_thresh;
	float kurt_thresh;
	float dm0_thresh;
	float cell_thresh;
	float gtest_thresh; // Threshold as far as the user is concerned - around 0.3
	float mean_max;
	float mean_min;
	float std_max;
	float std_min;
	int flag_grow;
	bool invert_freq;
	bool subtract_dm0;
	uint64_t interval_samps;
	int nf;
	int nt;
	int nbeams_per_ant; // number of beams*number of polarisations per antenna
	int nants; // number of antennas
	int nbits;
	int nsamps_per_int;
	bool polsum;
	DataOrder in_order;
	float global_scale;
};

class Rescaler {

public:

	RescaleOptions options;
	FreddaParams& params;
	RescaleOptions noflag_options;
	uint64_t num_flagged_times = 0;
	uint64_t num_flagged_beam_chans = 0;



	// Parameters
	/* I put so much effort into this I'm scared of deleting it now
	Rescaler(int _nbeams, int _nf, int _nt,
			float _target_mean, float _target_stdev, float _decay_constant,
			float _mean_thresh, float _std_thresh, float _kurt_thresh,
			float _dm0_thresh, float _cell_thresh,
			int _flag_grow, bool _invert_freq, bool _subtract_dm0);
			*/
	Rescaler(RescaleOptions& _options, FreddaParams& params);

	virtual ~Rescaler();

	void reset_output(array4d_t& rescale_buf); // Set output buffer to zero to start accumulating again

	void process_ant_block(void* read_buf, int iant, cudaStream_t stream=0);

	void finish_all_ants(array4d_t& outbuf); // Needs to be called after a cudaDeviceSynchronize

	void set_scaleoffset(float s_scale, float s_offset);

	void flag_channel(int channel); // Set weights to zero for all beams/antennas fo rthis channel
	int flag_frequencies_from_file(const char* filename); // Flag all frequencies in given file
	bool flag_frequency(float freq); // flag the channel with frequency nearst the given frequency. Ignored it out of band
	void flag_beam(int beamno);


private:

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
	array4d_t weights; // weights for antennas/beams.
	uint64_t sampnum;

	size_t in_buffer_bytes_per_ant;
	uint8_t* in_buffer_device;

	int num_elements;
	int num_elements_per_ant;

	CudaTimer tdump;
	int blocknum = 0;

	std::vector<Array4dDumper* > rescale_dumpers; // Dumpers that get updated per rescale interval
	std::vector<Array4dDumper* > block_dumpers; // Dumpers taht get updated per block

	// Need to separate this out for the first integration by antenna and don't write stats
	void update_rescale_parameters(RescaleOptions& options, int iant, cudaStream_t stream = 0);

	// Call after updatescaleoffset in the first run to reset to the beginning
	void reset_block_stats_for_first_block(int iant);

	// Update flagging statistics based on current block
	void update_flagging_statisics();
	void dump_rescale_data(); // Dump rescaler data to disk
	void dump_block_data(); // Dump block data to disk

	void apply_flags_and_sum(array4d_t& rescale_buf, RescaleOptions& options, int iant, cudaStream_t stream=0);

	template <int nsamps_per_word, typename wordT>
	void do_apply_flags_and_sum(array4d_t& rescale_buf, wordT* read_buf_device, RescaleOptions& options, int iant, cudaStream_t stream);


	void calculate_block_stats(void* read_buf_device, RescaleOptions &options, int iant, cudaStream_t stream=0);

	template <int nsamps_per_word, typename wordT>
	void do_calculate_block_stats(wordT* read_buf_device, RescaleOptions &options, int iant, cudaStream_t stream=0);


};

template <int nsamps_per_word, typename wordT>
void Rescaler::do_apply_flags_and_sum(array4d_t& rescale_buf, wordT* read_buf_device, RescaleOptions& options, int iant, cudaStream_t stream)
{
	int nbeams_in = options.nbeams_per_ant;
	int nf = rescale_buf.nx;
	int nt = rescale_buf.nz;
	int nwords = nf / nsamps_per_word;
	assert(nf % nsamps_per_word == 0);
	int boff = iant*options.nbeams_per_ant;
	assert(options.in_order == DataOrder::TFBP || options.in_order == DataOrder::BPTF);

	rescale_calc_dm0_kernel< nsamps_per_word, wordT > <<<nbeams_in, 256, 0, stream>>>(
			read_buf_device,
			offset.d_device,
			scale.d_device,
			dm0.d_device,
			dm0count.d_device,
			nf, nt,
			options.cell_thresh,
			options.in_order,
			boff);
	gpuErrchk(cudaPeekAtLastError());


	// Take the mean all the dm0 times into one big number per beam - this is the how we flag
	// short dropouts see ACES-209
	// probably could do this in rescale_calc_dm0_kernel after yu've done it
	// But i Haven't got htere yet.
	rescale_calc_dm0stats_kernel<<<1, nbeams_in, 0, stream>>>(
			dm0.d_device,
			dm0count.d_device,
			dm0stats.d_device,
			nt,
			boff);
	gpuErrchk(cudaPeekAtLastError());


	int nthread;
	if (nwords % (128/nsamps_per_word) == 0) {
		nthread = 128/nsamps_per_word; // Tuneable  -- needed if nf > 1024 which is a commmon limitiation for the number of threads.
	} else {
		nthread = nwords;
	}

	// TODO: CHECK blockdim doesn't exceed GPU capability
	dim3 griddim(nbeams_in, nwords/nthread);
	dim3 blockdim(nsamps_per_word, nthread);

	assert(griddim.y * blockdim.y == nwords);
	assert(blockdim.x == nsamps_per_word);
	assert(griddim.x == nbeams_in);

	rescale_apply_flags_and_add< nsamps_per_word, wordT ><<<griddim, blockdim, 0, stream>>>(
			read_buf_device,
			sum.d_device,
			sum2.d_device,
			sum3.d_device,
			sum4.d_device,
			decay_offset.d_device,
			nsamps.d_device,
			offset.d_device,
			scale.d_device,
			dm0.d_device,
			dm0count.d_device,
			dm0stats.d_device,
			rescale_buf.d_device,
			options.decay_constant,
			options.dm0_thresh,
			options.cell_thresh*options.target_stdev,
			nt,
			options.invert_freq,
			options.subtract_dm0,
			options.polsum,
			options.in_order,
			boff);

	gpuErrchk(cudaPeekAtLastError());

	sampnum += nt;

}


template <int nsamps_per_word, typename wordT>
void Rescaler::do_calculate_block_stats(wordT* read_buf_device, RescaleOptions &options, int iant, cudaStream_t stream)
{
	int nbeams_in = options.nbeams_per_ant;
	int nf = params.nf;
	int nt = params.nt;
	int nwords = nf / nsamps_per_word;
	assert(nf % nsamps_per_word == 0);
	int boff = iant*options.nbeams_per_ant;
	assert(options.in_order == DataOrder::TFBP || options.in_order == DataOrder::BPTF);

	rescale_calc_dm0_kernel< nsamps_per_word, wordT > <<<nbeams_in, 256, 0, stream>>>(
			read_buf_device,
			offset.d_device,
			scale.d_device,
			dm0.d_device,
			dm0count.d_device,
			nf, nt,
			options.cell_thresh,
			options.in_order,
			boff);
	gpuErrchk(cudaPeekAtLastError());


	// Take the mean all the dm0 times into one big number per beam - this is the how we flag
	// short dropouts see ACES-209
	// probably could do this in rescale_calc_dm0_kernel after yu've done it
	// But i Haven't got htere yet.
	rescale_calc_dm0stats_kernel<<<1, nbeams_in, 0, stream>>>(
			dm0.d_device,
			dm0count.d_device,
			dm0stats.d_device,
			nt,
			boff);
	gpuErrchk(cudaPeekAtLastError());


	int nthread;
	if (nwords % (128/nsamps_per_word) == 0) {
		nthread = 128/nsamps_per_word; // Tuneable  -- needed if nf > 1024 which is a commmon limitiation for the number of threads.
	} else {
		nthread = nwords;
	}

	// TODO: CHECK blockdim doesn't exceed GPU capability
	dim3 griddim(nbeams_in, nwords/nthread);
	dim3 blockdim(nsamps_per_word, nthread);

	assert(griddim.y * blockdim.y == nwords);
	assert(blockdim.x == nsamps_per_word);
	assert(griddim.x == nbeams_in);

	rescale_calc_stats< nsamps_per_word, wordT ><<<griddim, blockdim, 0, stream>>>(
			read_buf_device,
			sum.d_device,
			sum2.d_device,
			sum3.d_device,
			sum4.d_device,
			decay_offset.d_device,
			nsamps.d_device,
			offset.d_device,
			scale.d_device,
			dm0.d_device,
			dm0count.d_device,
			dm0stats.d_device,
			NULL, //rescale_buf.d_device,
			options.decay_constant,
			options.dm0_thresh,
			options.cell_thresh*options.target_stdev,
			nt,
			options.invert_freq,
			options.subtract_dm0,
			options.polsum,
			options.in_order,
			boff);

	gpuErrchk(cudaPeekAtLastError());
}




#endif /* RESCALER_H_ */
