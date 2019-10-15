#include "prepare.h"

/*
  Prepare Kernel Implementation 
  Arguments:
  in   (input)      --> Input RAW data
  cal   (input)     --> Input calibration
  sky   (input)     --> Input sky model
  out   (output)    --> Output RAW data
  average  (output) --> Output of average 
*/

// Order is assumed to be TBFP, BFP or BF
extern "C" {
  void prepare_knl(
		   const burst_datatype *in,  // Read-Only input raw data
		   const burst_datatype *cal, // Read-Only input calibration
		   const burst_datatype *sky, // Read-Only input sky model
		   burst_datatype *out,       // Output raw data
		   burst_datatype *average   // Output average data, lives with multiple kernel instances
		   )
  {
#pragma HLS INTERFACE m_axi port=in offset=slave bundle=gmem0
#pragma HLS INTERFACE m_axi port=cal offset=slave bundle=gmem1
#pragma HLS INTERFACE m_axi port=sky offset=slave bundle=gmem2
#pragma HLS INTERFACE m_axi port=out offset=slave bundle=gmem3
#pragma HLS INTERFACE m_axi port=average offset=slave bundle=gmem4
#pragma HLS INTERFACE s_axilite port=in  bundle=control
#pragma HLS INTERFACE s_axilite port=cal  bundle=control
#pragma HLS INTERFACE s_axilite port=sky bundle=control
#pragma HLS INTERFACE s_axilite port=out bundle=control
#pragma HLS INTERFACE s_axilite port=average bundle=control
#pragma HLS INTERFACE s_axilite port=return bundle=control

#pragma HLS DATA_PACK variable = in //struct_level
#pragma HLS DATA_PACK variable = cal //struct_level
#pragma HLS DATA_PACK variable = sky //struct_level
#pragma HLS DATA_PACK variable = out //struct_level
#pragma HLS DATA_PACK variable = average //struct_level

    int i;
    int j;
    int k;
    int l;
    int loc1;
    int loc2;
    
    burst_datatype in_temp[2];
    burst_datatype cal_temp[2];
    burst_datatype average_temp[2];
    burst_datatype sky_temp;
    burst_datatype out_temp;
    
#pragma HLS DATA_PACK variable = in_temp
#pragma HLS DATA_PACK variable = cal_temp
#pragma HLS DATA_PACK variable = average_temp
#pragma HLS DATA_PACK variable = sky_temp
#pragma HLS DATA_PACK variable = out_temp
    
#pragma HLS array_partition variable=in_temp complete
#pragma HLS array_partition variable=cal_temp complete
#pragma HLS array_partition variable=average_temp complete
#pragma HLS dataflow
    
    for(i = 0; i < NBASELINE; i++){
      for(j = 0; j < NBURST; j++){
      cal_loc1: loc1 = i*NBURST + j;
	
      read_cal0: cal_temp[0] = cal[2*loc1];
      read_cal1: cal_temp[1] = cal[2*loc1+1];	
      read_sky: sky_temp     = sky[loc1];
	
      reset_average: for(l = 0; l < BURST_DATA_SIZE; l++){
	  average_temp[0].data[l].real(0);
	  average_temp[0].data[l].imag(0);
	  average_temp[1].data[l].real(0);
	  average_temp[1].data[l].imag(0);
	}
	
	for(k = 0; k < NTIME_PER_BUFBLOCK; k++){
#pragma HLS PIPELINE II=1
	cal_loc2: loc2 = k*NBASELINE*NBURST + i*NBURST + j;
	  
	read_in0: in_temp[0] = in[2*loc2];
	read_in1: in_temp[1] = in[2*loc2+1];
	  
	  for(l = 0; l < BURST_DATA_SIZE; l++){
	  do_average0: average_temp[0].data[l]+=in_temp[0].data[l];
	  do_average1: average_temp[1].data[l]+=in_temp[1].data[l];
	  }
	  
	  for(l = 0; l < BURST_DATA_HALF; l++){
	  do_out0: out_temp.data[l] = in_temp[0].data[2*l] * cal_temp[0].data[2*l] +
	      in_temp[0].data[2*l+1] * cal_temp[0].data[2*l+1] -
	      sky_temp.data[l];
	    
	  do_out1:out_temp.data[BURST_DATA_HALF+l] = in_temp[1].data[2*l] * cal_temp[1].data[2*l] +
	      in_temp[1].data[2*l+1] * cal_temp[1].data[2*l+1] -
	      sky_temp.data[BURST_DATA_HALF+l];	    
	  }
	write_out: out[loc2] = out_temp;
	}
      write_average0: average[2*loc1] = average_temp[0];	
      write_average1: average[2*loc1+1] = average_temp[1];
      }
      
    }
  }
}
