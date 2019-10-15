/*******************************************************************************
** C Code to represent kernel, not kernel code
*******************************************************************************/

#include "util_sdaccel.h"
#include "prepare.h"

int prepare(cl_complex *in, cl_complex *cal, cl_complex *sky, cl_complex *out, cl_complex *average)
{
  cl_uint i;
  cl_uint j;
  cl_uint k;
  cl_uint loc1;
  cl_uint loc2;

  for(i=0; i<NBASELINE; i++){
    for(j=0; j<NCHAN; j++){
      loc1 = i * NCHAN + j;      
      for(k=0; k<NTIME_PER_BUFBLOCK; k++){
	
	loc2 = k * NCHAN * NBASELINE + i * NCHAN + j;      
	
	average[loc1*2]   += in[loc2*2];
	average[loc1*2+1] += in[loc2*2+1];
	
	out[loc2] = in[loc2*2] * cal[loc1*2] +
	  in[loc2*2+1] * cal[loc1*2+1] -
	  sky[loc1];
      }
    }
  }
  
  return EXIT_SUCCESS;
}
