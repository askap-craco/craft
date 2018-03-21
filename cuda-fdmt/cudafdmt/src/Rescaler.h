/*
 * Rescaler.h
 *
 *  Created on: 21 Mar 2018
 *      Author: ban115
 */

#ifndef RESCALER_H_
#define RESCALER_H_

#include "rescale.h"
#include "array.h"

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
	int flag_grow;
	bool invert_freq;
	bool subtract_dm0;
	uint64_t interval_samps;
	int nf;
	int nt;
	int nbeams;
	int nbits;
};

class Rescaler {

public:
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
	uint64_t sampnum;
	int num_elements;
	RescaleOptions options;


	// Parameters
	/* I put so much effort into this I'm scared of deleting it now
	Rescaler(int _nbeams, int _nf, int _nt,
			float _target_mean, float _target_stdev, float _decay_constant,
			float _mean_thresh, float _std_thresh, float _kurt_thresh,
			float _dm0_thresh, float _cell_thresh,
			int _flag_grow, bool _invert_freq, bool _subtract_dm0);
			*/
	Rescaler(RescaleOptions& _options);

	virtual ~Rescaler();

	void update_scaleoffset();
	void set_scaleoffset(float s_scale, float s_offset);
	void update_and_transpose(array4d_t& rescale_buf, void* read_buf_device);

	template <int nsamps_per_word, typename wordT>
	void do_update_and_transpose(array4d_t& rescale_buf, wordT* read_buf_device);
};

/*
template <int nsamps_per_word, class wordT>
Rescaler<nsamps_per_word, wordT>::Rescaler(int _nbeams, int _nf, int _nt,
		float _target_mean, float _target_stdev, float _decay_constant,
		float _mean_thresh, float _std_thresh, float _kurt_thresh,
		float _dm0_thresh, float _cell_thresh,
		int _flag_grow, bool _invert_freq, bool _subtract_dm0) : nbeams(_nbeams), nf(_nf), nt(_nt),
		target_mean(_target_mean),
		target_stdev(_target_stdev), decay_constant(_decay_constant),
		mean_thresh(_mean_thresh), std_thresh(_std_thresh), kurt_thresh(_kurt_thresh),
		dm0_thresh(_dm0_thresh), cell_thresh(_cell_thresh),
		flag_grow(_flag_grow), invert_freq(_invert_freq), subtract_dm0(_subtract_dm0) {
		*/



template <int nsamps_per_word, typename wordT>
void Rescaler::do_update_and_transpose(array4d_t& rescale_buf, wordT* read_buf_device)
{
	int nbeams = rescale_buf.nw;
	int nf = rescale_buf.nx;
	int nt = rescale_buf.nz;
	int nwords = nf / nsamps_per_word;
	assert(nf % nsamps_per_word == 0);

	// clear output
	array4d_cuda_memset(&rescale_buf, 0);

	rescale_calc_dm0_kernel< nsamps_per_word, wordT > <<<nbeams, 256>>>(
			read_buf_device,
			offset.d_device,
			scale.d_device,
			dm0.d_device,
			dm0count.d_device,
			nf, nt,
			options.cell_thresh);

	// Take the mean all the dm0 times into one big number per beam - this is the how we flag
	// short dropouts see ACES-209
	// probably could do this in rescale_calc_dm0_kernel after yu've done it
	// But i Haven't got htere yet.
	rescale_calc_dm0stats_kernel<<<1, nbeams>>>(
			dm0.d_device,
			dm0count.d_device,
			dm0stats.d_device,
			nt);

	dim3 blockdim(nsamps_per_word, nwords);

	rescale_update_and_transpose_float_kernel< nsamps_per_word, wordT ><<<nbeams, blockdim>>>(
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
			options.subtract_dm0);


	sampnum += nt;
	gpuErrchk(cudaDeviceSynchronize());
}




#endif /* RESCALER_H_ */
