#ifndef _RESCALE_H
#define _RESCALE_H

typedef struct {
  /* stuff for rescaling */
  float* sum;
  float* sumsq;
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

rescale_t* rescale_allocate(rescale_t* rescale, uint64_t nelements) ;

void rescale_update_scaleoffset(rescale_t* rescale);
void rescale_update_none(rescale_t* rescale, float* inx, float*outx);
void rescale_update_float(rescale_t* rescale, float* fdata, float* sampbuf, unsigned istart);
void rescale_update_float_polsum(rescale_t* rescale, float* fdata, float* sampbuf, unsigned istart);
void rescale_update_uint8(rescale_t* rescale, float* fdata, uint8_t* sampbuf, unsigned istart);
void rescale_update_uint8_polsum(rescale_t* rescale, float* fdata, uint8_t* sampbuf, unsigned istart);
void rescale_update_int8(rescale_t* rescale, float* in, int8_t*  out);
void rescale_update_decay_float(rescale_t* rescale, float* in, float* out);
void rescale_update_decay_uint8(rescale_t* rescale, float* in, uint8_t* out);

#endif
