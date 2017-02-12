#ifndef _RESCALE_H
#define _RESCALE_H

#include <stdint.h>
#include "array.h"


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
void rescale_update_and_transpose_float_gpu(rescale_gpu_t& rescale, array4d_t& rescale_buf, const uint8_t* read_buf, bool invert_freq);


#endif
