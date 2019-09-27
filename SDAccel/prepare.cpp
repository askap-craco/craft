/**********
Copyright (c) 2018, Xilinx, Inc.
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors
may be used to endorse or promote products derived from this software
without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********/

/*******************************************************************************
Description:
    HLS pragmas can be used to optimize the design : improve throughput, reduce latency and 
    device resource utilization of the resulting RTL code
    This is vector addition example to demonstrate how HLS optimizations are used in kernel. 
*******************************************************************************/


#include "prepare.h"
#include <complex>

/*
    Prepare Kernel Implementation 
    Arguments:
        in   (input)      --> Input RAW data
        cal   (input)     --> Input calibration
	sky   (input)     --> Input sky model
        out   (output)    --> Output RAW data
        average  (output) --> Output of average 
   */

//typedef prepare_t std::complex<float>;

// Order is assumed to be TBCP, BCP or BC
extern "C" {
void prepare(
        const float *in,  // Read-Only input raw data
        const float *cal, // Read-Only input calibration 
        const float *sky, // Read-Only input sky model
	float *out,       // Output raw data
	float *average   // Output average data, lives with multiple kernel instances
        )
{
  // SDAccel kernel must have one and only one s_axilite interface which will be used by host application to configure the kernel.
  // Here bundle control is defined which is s_axilite interface and associated with all the arguments (in1, in2, out and size),
  // control interface must also be associated with "return".
  // All the global memory access arguments must be associated to one m_axi(AXI Master Interface). Here all three arguments(in1, in2, out) are 
  // associated to bundle gmem which means that a AXI master interface named "gmem" will be created in Kernel and all these variables will be 
  // accessing global memory through this interface.
  // Multiple interfaces can also be created based on the requirements. For example when multiple memory accessing arguments need access to
  // global memory simultaneously, user can create multiple master interfaces and can connect to different arguments.
#pragma HLS INTERFACE m_axi port=in  offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=cal  offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=sky offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=out offset=slave bundle=gmem
#pragma HLS INTERFACE m_axi port=average offset=slave bundle=gmem
#pragma HLS INTERFACE s_axilite port=in  bundle=control
#pragma HLS INTERFACE s_axilite port=cal  bundle=control
#pragma HLS INTERFACE s_axilite port=sky bundle=control
#pragma HLS INTERFACE s_axilite port=out bundle=control
#pragma HLS INTERFACE s_axilite port=average bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control
  
  float cal_buf[NCHAN*4];
  float sky_buf[NCHAN*2];
  float in_buf[NCHAN*4];
  float out_buf[NCHAN*2];
  float average_buf[NCHAN*4]={0};
  int i; 
  int j;
  int k;
  int loc;

  for(i = 0; i < NBASELINE; i++)
    {
      //Burst read in calibration 
    read_cal: for(j = 0; j < NCHAN; j++) {
	loc = i * NCHAN + j;
	cal_buf[4*j]     = cal[4*loc];
	cal_buf[4*j + 1] = cal[4*loc + 1];
	cal_buf[4*j + 2] = cal[4*loc + 2];
	cal_buf[4*j + 3] = cal[4*loc + 3];
      }
      
      // Burst read in sky model
    read_sky: for(j = 0; j < NCHAN; j++) {
	loc = i * NCHAN + j;
	sky_buf[2*j]     = sky[2*loc];
	sky_buf[2*j + 1] = sky[2*loc + 1];
      }
      
    reset_average: for(j = 0; j < NCHAN; j++){
	average_buf[j*4]    =0;
	average_buf[j*4 + 1]=0;
	average_buf[j*4 + 2]=0;
	average_buf[j*4 + 3]=0;
      }
      
      for(j = 0; j < NTIME_PER_BUFBLOCK; j++) {
      read_in:for(k = 0; k < NCHAN; k++){
	  loc = j * NBASELINE * NCHAN + i * NCHAN +  k;
	  in_buf[k*4]     = in[loc*4];
	  in_buf[k*4 + 1] = in[loc*4 + 1];
	  in_buf[k*4 + 2] = in[loc*4 + 2];
	  in_buf[k*4 + 3] = in[loc*4 + 3];     
	}
	
      prepare: for(k = 0; k < NCHAN; k++){
	  // Do the calculation
	  average_buf[k*4]    +=in_buf[k*4];
	  average_buf[k*4 + 1]+=in_buf[k*4 + 1];
	  average_buf[k*4 + 2]+=in_buf[k*4 + 2];
	  average_buf[k*4 + 3]+=in_buf[k*4 + 3];
	  
	  out_buf[k*2] = in_buf[k*4]*cal_buf[k*4] - in_buf[k*4 + 1]*cal_buf[k*4 + 1] + // Real of pol 1 after cal
	    in_buf[k*4 + 2]*cal_buf[k*4 + 2] - in_buf[k*4 + 3]*cal_buf[k*4 + 3] - // Real of pol 2 after cal
	    sky_buf[k*2];// Real of sky
	  
	  out_buf[k*2 + 1] = in_buf[k*4]*cal_buf[k*4 + 1] + in_buf[k*4 + 1]*cal_buf[k*4] + // Image of pol 1 after cal
	    in_buf[k*4 + 2]*cal_buf[k*4 + 3] + in_buf[k*4 + 3]*cal_buf[k*4 + 2] - // Image of pol 2 after cal
	    sky_buf[k*2 + 1];// Image of sky	  
	}
    
	// Burst output raw data
      write_out: for(k = 0; k < NCHAN; k++){
	  loc = j * NBASELINE * NCHAN + i * NCHAN +  k;
	  out[loc*2]     = out_buf[k*2];
	  out[loc*2 + 1] = out_buf[k*2 + 1];
	}  
      }      
      // Burst output average data
    write_average: for(j = 0; j < NCHAN; j++) {
	loc = i * NCHAN + j;
	average[4*loc]     = average_buf[4*j];
	average[4*loc + 1] = average_buf[4*j + 1];
	average[4*loc + 2] = average_buf[4*j + 2];
	average[4*loc + 3] = average_buf[4*j + 3];
      }
    }
}
}
