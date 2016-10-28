/*
 * Rescaling utilities
 * Author: Keith Bannister <keith.bannister@csiro.au>
 */

#include "rescale.h"
#include <stdint.h>
#include <unistd.h>
#include <assert.h>
#include <math.h>

void* rescale_malloc(size_t sz)
{
	void* ptr = malloc(sz);
	assert(ptr);
	return ptr;
}

void rescale_update_scaleoffset(rescale_t* rescale)
{
  assert(rescale->interval_samps >= 0);
  assert(rescale->target_stdev > 0);
  float nsamp = (float) rescale->sampnum;
  for (unsigned i = 0; i < rescale->num_elements; i++) {
    float mean = rescale->sum[i]/nsamp;
    float meansq = rescale->sumsq[i]/nsamp;
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
    rescale->sumsq[i] = 0.0;
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
    rescale->sumsq[i] += vin2;
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
    rescale->sumsq[j] += vin*vin;
    
    float uin = inx[j+1];
    rescale->sum[j+1] += uin;
    rescale->sumsq[j+1] += uin*uin;

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
    rescale->sumsq[i] += vin*vin;
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
    rescale->sumsq[j] += vin*vin;

    float uin=in[j+1];
    rescale->sum[j+1] += uin;
    rescale->sumsq[j+1] += uin*uin;

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
    rescale->sumsq[i] += vin*vin;
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
 
void rescale_update_decay_float(rescale_t* rescale, float* __restrict__ in, float* __restrict__ out)
{
  float k = rescale->decay_constant;

  for (unsigned i = 0; i < rescale->num_elements; i++) {
    float vin = in[i];
    rescale->sum[i] += vin;
    rescale->sumsq[i] += vin*vin;
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
    rescale->sumsq[i] += vin*vin;
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
  rescale->sumsq = (float*) rescale_malloc(sz);
  rescale->scale = (float*) rescale_malloc(sz);
  rescale->offset = (float*) rescale_malloc(sz);
  rescale->decay_offset = (float*) rescale_malloc(sz);
  rescale->sampnum = 0;
  rescale->num_elements = nelements;
  
  //rescale_update_scaleoffset(rescale);
  return rescale;

}
