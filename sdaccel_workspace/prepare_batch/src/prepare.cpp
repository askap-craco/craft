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
	#pragma HLS INTERFACE s_axilite port=in bundle=control
	#pragma HLS INTERFACE s_axilite port=cal bundle=control
	#pragma HLS INTERFACE s_axilite port=sky bundle=control
	#pragma HLS INTERFACE s_axilite port=out bundle=control
	#pragma HLS INTERFACE s_axilite port=average bundle=control
	#pragma HLS INTERFACE s_axilite port=return bundle=control

	#pragma HLS DATA_PACK variable = in
	#pragma HLS DATA_PACK variable = cal
	#pragma HLS DATA_PACK variable = sky
	#pragma HLS DATA_PACK variable = out
	#pragma HLS DATA_PACK variable = average

	int i;
    int j;
    int k;
    int l;
    int loc;
    burst_datatype cal_temp[2*NBURST];
    burst_datatype in_temp[2*NBURST];
    burst_datatype sky_temp[NBURST];
    burst_datatype out_temp[NBURST];
    burst_datatype average_temp[2*NBURST];

    #pragma HLS RESOURCE variable=cal_temp core=XPM_MEMORY uram // utilizes URAM
	#pragma HLS RESOURCE variable=in_temp core=XPM_MEMORY uram // utilizes URAM
	#pragma HLS RESOURCE variable=sky_temp core=XPM_MEMORY uram // utilizes URAM
	#pragma HLS RESOURCE variable=out_temp core=XPM_MEMORY uram // utilizes URAM
    #pragma HLS RESOURCE variable=average_temp core=XPM_MEMORY uram // utilizes URAM

    /*
    #pragma HLS array_partition variable=in_temp complete
	#pragma HLS array_partition variable=cal_temp complete
	#pragma HLS array_partition variable=average_temp complete
	#pragma HLS array_partition variable=out_temp complete
	#pragma HLS array_partition variable=sky_temp complete
    */
	#pragma HLS dataflow

    for(i = 0; i < NBASELINE; i++){

    	read_cal: for(j = 0; j < 2*NBURST; j++){
    		#pragma HLS PIPELINE II=1
    		loc = 2*i*NBURST + j;
    		cal_temp[j] = cal[loc];
    	}

    	read_sky: for(j = 0; j < NBURST; j++){
			#pragma HLS PIPELINE II=1
    		loc = i*NBURST + j;
    		sky_temp[j] = sky[loc];
    	}

    	reset_average: for(j = 0; j < 2*NBURST; j++){
    		for(l = 0; l < BURST_DATA_SIZE; l++){
    			#pragma HLS PIPELINE II=1
    			average_temp[j].data[l].real(0);
    			average_temp[j].data[l].imag(0);
    		}
    	}

    for(k = 0; k < NTIME_PER_BUFBLOCK; k++){
    	read_in: for(j = 0; j < 2*NBURST; j++){
    		#pragma HLS PIPELINE II=1
    		loc = 2*k*NBASELINE*NBURST + 2*i*NBURST + j;
    		in_temp[j] = in[loc];
    	}

    	do_average: for(j = 0; j < 2*NBURST; j++){
    		for(l = 0; l < BURST_DATA_SIZE; l++){
    			#pragma HLS PIPELINE II=1
    			average_temp[j].data[l] += in_temp[j].data[l];
    		}
    	}
	
    	do_out: for(j = 0; j < NBURST; j++){
    		for(l = 0; l < BURST_DATA_HALF; l++){
    			#pragma HLS PIPELINE II=1
    			do_out0: out_temp[j].data[l] = in_temp[2*j].data[2*l] * cal_temp[2*j].data[2*l] +
    					in_temp[2*j].data[2*l+1] * cal_temp[2*j].data[2*l+1] -
						sky_temp[j].data[l];
	    
    			do_out1:out_temp[j].data[BURST_DATA_HALF+l] = in_temp[2*j+1].data[2*l] * cal_temp[2*j+1].data[2*l] +
    					in_temp[2*j+1].data[2*l+1] * cal_temp[2*j+1].data[2*l+1] -
						sky_temp[j].data[BURST_DATA_HALF+l];
    		}
    	}

    	write_out: for(j = 0; j < NBURST; j++) {
    		#pragma HLS PIPELINE II=1
    		loc = k*NBASELINE*NBURST + i*NBURST + j;
    		out[loc] = out_temp[j];
    	}
    }

    write_average: for (j = 0; j < 2*NBURST; j++){
    	#pragma HLS PIPELINE II=1
    	loc = 2*i*NBURST + j;
    	average[loc] = average_temp[j];
    }
    }
}
}
