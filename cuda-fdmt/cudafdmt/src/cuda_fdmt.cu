/*
 * cuda_fdmt.cu
 *
 *  Created on: 4 Oct 2016
 *      Author: ban115
 */
#include <assert.h>
#import "fdmt.h"
#import "cuda_fdmt.h"

//#import "gpu_kernels.h"

void __global__ cuda_fdmt_iteration1(fdmt_t fdmt,
		int iteration_num,
		array4d_t indata,
		array4d_t outdata)
{
	float df = fdmt.df; // channel resolution
	float delta_f = (float)(1 << iteration_num) * df; // Resolution of current iteration
	int delta_t = calc_delta_t(&fdmt, fdmt.fmin, fdmt.fmin+delta_f); // Max DM

	// Outdata has size (nbeams, o_nf, o_nd1, fdmt.nt)
	int onx = indata.nx/2 + indata.nx % 2;
	outdata.nw = indata.nw;
	outdata.nx = indata.nx/2 + indata.nx % 2; // Add 1 to the frequency dimension if it's not divisible by 2
	outdata.ny = delta_t + 1;
	outdata.nz = indata.nz;


	//
	//	assert(array4d_size(&outdata) <= fdmt.state_size);
	//
	//	//    printf("iteration %d df %f delta_f %f delta_t %d output nx=%d ny=%d nz%d\n",
	//	//           iteration_num, df, delta_f, delta_t, outdata.nx, outdata.ny, outdata.nz);
	//
	//	// zero that output baby
	//	//bzero(outdata.d, outdata.nw*outdata.nx * outdata.ny * outdata.nz * sizeof(fdmt_dtype));
	//	//array4d_cuda_memset(&outdata, 0);
	//
	int shift_input = 0; // ?
	int shift_output = 0; // ?

	float fjumps = (float)outdata.nx; // Output number of channels
	float frange = fdmt.fmax - fdmt.fmin; // Width of band
	float fmin = fdmt.fmin; // Bottom of band

	float correction = 0.0;
	if (iteration_num > 0) {
		correction = df/2.0;
	}
	//
	//	//assert(indata.nw == fdmt.nbeams);
	////	// For each output sub-band
	for (int iif = 0; iif < onx; iif++) {
		float f_start = frange/fjumps * (float)iif + fmin; // Top freq of subband
		float f_end = frange/fjumps*((float)iif + 1) + fmin; // Bottom freq of subband
		float f_middle = (f_end - f_start)/2.0 + f_start - correction; // Middle freq of subband, less 0.5xresolution
		float f_middle_larger = (f_end - f_start)/2.0 + f_start + correction; // Middle freq of subband + 0.5x resolution (helps with rounding)
		//
		//		// Max DM for this subband
		//		int delta_t_local = calc_delta_t(&fdmt, f_start, f_end) + 1;

		//		// For each DM relevant for this subband
		//		for (int idt = 0; idt < delta_t_local; idt++) {
		//			int dt_middle = roundf(idt * cff(f_middle, f_start, f_end, f_start)); // Dt for middle freq less 0.5xresolution
		//			int dt_middle_index = dt_middle + shift_input;
		//			int dt_middle_larger = roundf(idt * cff(f_middle_larger, f_start, f_end, f_start)); // Dt for middle freq +0.5x resolution
		//			int dt_rest = idt - dt_middle_larger;
		//			int dt_rest_index = dt_rest + shift_input;
		//
		//			int itmin = 0;
		//			int itmax = dt_middle_larger;
		//
		//			//Output[i_F,i_dT + ShiftOutput,i_T_min:i_T_max] = Input[2*i_F, dT_middle_index,i_T_min:i_T_max];
		//			//int outidx = array4d_idx(outdata, beam, iif, idt+shift_output, 0);
		//			//int inidx1  = array4d_idx(indata, beam, 2*iif, dt_middle_index, 0);
		//
		//
		//			//printf("iteration %d channel %d freq %f idt %d dt_local "
		//			//	  "%d dt_middle %d dt_middle_larger %d dt_rest %d\n",
		//			//	  iteration_num, iif, f_middle, idt, delta_t_local, dt_middle_index, dt_middle_larger, dt_rest_index);
		//
		//			// Here we handle the edge effects and set
		//			// OUtput state[freq, idx, 0:dtmin] = input_state[2xfreq, dt_middle, 0:dtmin]
		//			// where the DM would have overun the available times
		//			// This needs to be fixed for more careful time overlapping
		//			coord3_t dst_start = {.x = iif, .y = idt+shift_output, .z = 0};
		//			coord3_t src1_start = {.x = 2*iif, .y = dt_middle_index, .z = 0};
		//			//array_gpu_copy1(&outdata, &indata, &dst_start, &src1_start, dt_middle_larger);
		//			//cpu_copy2(&outdata.d[outidx + itmin], &indata.d[inidx1 + itmin], (itmax - itmin));
		//			//for (int i = itmin; i < itmax; i++) {
		//			//outdata.d[outidx + i] = indata.d[inidx1 + i];
		//			//}
		//
		//
		//			// Now we work on the remaining times that are guaranteed not to overrun the input dimensions
		//			itmin = dt_middle_larger;
		//			itmax = fdmt.max_dt;
		//
		//			coord3_t src2_start = {.x = 2*iif + 1, .y = dt_rest_index, .z = 0};
		//			// src and dst now start from a bit offset
		//			src1_start.z = dt_middle_larger;
		//			dst_start.z = dt_middle_larger;
		//			int zcount = itmax - itmin;
		//
		//
		//			if (2*iif + 1 < indata.nx) { // If the input data has this channel, we'll add it in
		//				//Output[i_F,i_dT + ShiftOutput,i_T_min:i_T_max] = Input[2*i_F, dT_middle_index,i_T_min:i_T_max] + Input[2*i_F+1, dT_rest_index,i_T_min - dT_middle_larger:i_T_max-dT_middle_larger]
		//				// playinga trick here - we're always addign the fastest moving index
		//				// Putting -dt_middle_larger in array3d_idx would have caused an assertion failure
		//				// But ofsetting by dt_middle_larger at the end, we get the best of all worlds
		//
		//				//int inidx2 = array4d_idx(indata, beam, 2*iif+1, dt_rest_index, 0) - dt_middle_larger;
		//
		//				//array_gpu_sum1(&outdata, &indata, &dst_start, &src1_start, &src2_start, zcount);
		//
		//				//for(int i = itmin; i < itmax; i++) {
		//				//  outdata.d[outidx + i] = indata.d[inidx1 + i] + indata.d[inidx2 + i];
		//				//}
		//				//cpu_sum1(&outdata.d[outidx + itmin], &indata.d[inidx1+itmin], &indata.d[inidx2+itmin], itmax-itmin);
		//
		//
		//			} else { // Just copy the input over. which basically assumes the upper channel is flaggedd/0
		//				// TODO: Could probably be done outside the iif loop to save evalutating IFs, but
		//				// Too tricky for the moment.
		//				//cpu_copy2(&outdata.d[outidx + itmin], &indata.d[inidx1 + itmin], (itmax - itmin));
		//				//for(int i = itmin; i < itmax; i++) {
		//				//  outdata.d[outidx + i] = indata.d[inidx1 + i];
		//				//	}
		//				//array_gpu_copy1(&outdata, &indata, &dst_start, &src1_start, zcount);
		//			}
		//		}

	}

}


