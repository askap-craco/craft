/*******************************************************************************
** C Code to represent kernel, not kernel code
*******************************************************************************/

#include "util_sdaccel.h"
#include "prepare.h"

int prepare(float *in, float *cal, float *sky, float *out, float *average)
{
  cl_uint i;
  cl_uint j;
  cl_uint k;
  cl_uint loc_raw;
  cl_uint loc_average;
  
  for(i=0; i<NTIME_PER_BUFBLOCK; i++){
    for(j=0; j<NCHAN; j++){
      for(k=0; k<NBASELINE; k++){
	loc_raw = i * NCHAN * NBASELINE + j * NBASELINE + k;
	loc_average = j * NBASELINE + k;

	/* Get average data */
	average[loc_average*4]+=in[loc_raw*4];
	average[loc_average*4 + 1]+=in[loc_raw*4 + 1];
	average[loc_average*4 + 2]+=in[loc_raw*4 + 2];
	average[loc_average*4 + 3]+=in[loc_raw*4 + 3];
	
	/* Get output raw data */
	out[loc_raw*2] = in[loc_raw*4]*cal[loc_average*4] - in[loc_raw*4 + 1]*cal[loc_average*4 + 1] + // Real of pol 1 after cal
	  in[loc_raw*4 + 2]*cal[loc_average*4 + 2] - in[loc_raw*4 + 3]*cal[loc_average*4 + 3] - // Real of pol 2 after cal
	  sky[loc_average*2];// Real of sky
	
	out[loc_raw*2 + 1] = in[loc_raw*4]*cal[loc_average*4 + 1] + in[loc_raw*4 + 1]*cal[loc_average*4] + // Image of pol 1 after cal
	  in[loc_raw*4 + 2]*cal[loc_average*4 + 3] + in[loc_raw*4 + 3]*cal[loc_average*4 + 2] - // Image of pol 2 after cal
	  sky[loc_average*2 + 1];// Image of sky	  
      }
    }
  }

  
  return EXIT_SUCCESS;
}
