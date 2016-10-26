#include <iostream>
#include "cpu_kernels.h"
#include "gpu_kernels.h"
//#include "cuda_fdmt.h"
#include "fdmt.h"
#include "CudaTimer.h"


float dm_delay(const float f1, const float f2) {
	return 4.14e9f*(isquaref(f1) - isquaref(f2));
}

__host__ __device__ float squaref(const float f)
{
	return f*f;
}

__host__ __device__ float isquaref(const float f)
{
	return 1.0f/(f*f);
}

__host__ __device__ float cff(float f1_start, float f1_end, float f2_start, float f2_end)
{
	float rf = (isquaref(f1_start) - isquaref(f1_end))/(isquaref(f2_start) - isquaref(f2_end));

	return rf;
}

__device__ inline float cff2(const float f1_start, const float f1_end, const float f2_start)
{
	//	float if1 = __frcp_rz(f1_start*f1_start);
	//	float if2 = __frcp_rz(f1_end*f1_end);
	//	float if3 = __frcp_rz(f2_start*f2_start);
	//
	//	float rf = (if1 - if2)/(if3 - if2);

	return 0.0f;
}

__host__ __device__ int calc_delta_t(const fdmt_t* fdmt, float f_start, float f_end)
{
	float rf = cff(f_start, f_end, fdmt->fmin, fdmt->fmax);
	float delta_tf = ((float)fdmt->max_dt-1.0) * rf;
	int delta_t = (int)ceilf(delta_tf);

	//  printf("delta t: rf %f delta_tf %f delta_t %d\n", rf, delta_tf, delta_t);
	return delta_t;
}

__host__ FdmtIteration* fdmt_save_iteration(fdmt_t* fdmt, const int iteration_num, const array4d_t* indata, array4d_t* outdata)
{
	float df = fdmt->df; // channel resolution
	float delta_f = (float)(1 << iteration_num) * df; // Resolution of current iteration
	int delta_t = calc_delta_t(fdmt, fdmt->fmin, fdmt->fmin+delta_f); // Max DM

	//FdmtIteration* iter = &fdmt->iterations[iteration_num-1];
	int nf = indata->nx/2 + indata->nx % 2; // Add 1 to the frequency dimension if it's not divisible by 2
	int ndt = delta_t + 1;

	// Outdata has size (nbeams, o_nf, o_nd1, fdmt->nt)
	outdata->nw = indata->nw;
	outdata->nx = nf; // Add 1 to the frequency dimension if it's not divisible by 2
	outdata->ny = ndt;
	outdata->nz = indata->nz;

	assert(array4d_size(outdata) <= fdmt->state_size);

	FdmtIteration* iter = new FdmtIteration(outdata->nw, outdata->nx, outdata->ny, outdata->nz);
	fdmt->iterations.push_back(iter);
	float fjumps = (float)nf; // Output number of channels
	float frange = fdmt->fmax - fdmt->fmin; // Width of band
	float fmin = fdmt->fmin; // Bottom of band

	float correction = 0.0;
	if (iteration_num > 0) {
		correction = df/2.0;
	}

	int shift_input = 0;
	int shift_output = 0;


	// For each output sub-band
	for (int iif = 0; iif < nf; iif++) {
		float f_start = frange/fjumps * (float)iif + fmin; // Top freq of subband
		float f_end = frange/fjumps*((float)iif + 1) + fmin; // Bottom freq of subband
		float f_middle = (f_end - f_start)/2.0 + f_start - correction; // Middle freq of subband, less 0.5xresolution
		float f_middle_larger = (f_end - f_start)/2.0 + f_start + correction; // Middle freq of subband + 0.5x resolution (helps with rounding)

		// Max DM for this subband
		int delta_t_local = calc_delta_t(fdmt, f_start, f_end) + 1;

		iter->add_subband(delta_t_local);

		// For each DM relevant for this subband
		for (int idt = 0; idt < delta_t_local; idt++) {
			// calculate parameters that take ages to do on the GPU for reasons not yet understood
			int dt_middle = roundf(idt * cff(f_middle, f_start, f_end, f_start)); // Dt for middle freq less 0.5xresolution
			int dt_middle_index = dt_middle + shift_input;
			int dt_middle_larger = roundf(idt * cff(f_middle_larger, f_start, f_end, f_start)); // Dt for middle freq +0.5x resolution
			int dt_rest = idt - dt_middle_larger;
			int dt_rest_index = dt_rest + shift_input;

			int itmin = 0;
			int itmax = dt_middle_larger;

			coord3_t dst_start = {.x = iif, .y = idt+shift_output, .z = 0};
			coord3_t src1_start = {.x = 2*iif, .y = dt_middle_index, .z = 0};
			coord3_t src2_start = {.x = 2*iif + 1, .y = dt_rest_index, .z = 0};

			//array_gpu_copy1(outdata, indata, &dst_start, &src1_start, dt_middle_larger);

			// Now we work on the remaining times that are guaranteed not to overrun the input dimensions
			itmin = dt_middle_larger;
			itmax = fdmt->max_dt;

			// src and dst now start from a bit offset
			src1_start.z = dt_middle_larger;
			dst_start.z = dt_middle_larger;
			int zcount = itmax - itmin;

			int mint = dt_middle_larger;
			int src1_offset = array4d_idx(indata, 0, 2*iif, dt_middle_index, 0);
			int src2_offset = array4d_idx(indata, 0, 2*iif+1, dt_rest_index, 0) - mint;
			int out_offset = array4d_idx(outdata, 0, iif, idt, 0);
			//			printf("iter %d iif %03d idt %02d src1_off %06d src2_off %06d out_off %06d maxt %02d dtmid %d dtr %d dtmidlg %d in [%d,%d,%d,%d]\n",
			//					iteration_num, iif, idt, src1_offset, src2_offset, out_offset, maxt, dt_middle_index, dt_rest_index, dt_middle_larger, indata->nw, indata->nx, indata->ny, indata->nz);

			iter->save_subband_values(idt, src1_offset, src2_offset, out_offset, mint);
		}
	}
	iter->copy_to_device();

	return iter;
}


int fdmt_create(fdmt_t* fdmt, float fmin, float fmax, int nf, int max_dt, int nt, int nbeams)
{
	fdmt->max_dt = max_dt;
	fdmt->nt = nt;
	fdmt->fmin = fmin;
	fdmt->fmax = fmax;
	fdmt->nf = nf;
	fdmt->df = (fdmt->fmax - fdmt->fmin)/((float) fdmt->nf);
	fdmt->order = (int)ceil(log(fdmt->nf)/log(2.0));
	fdmt->nbeams = nbeams;
	assert(nf > 0);
	assert(max_dt > 0);
	assert(nt > 0);
	assert(fdmt->max_dt >= fdmt->nt);
	assert(fdmt->max_dt % fdmt->nt == 0); // max_dt needs to be a multipel of nt
	assert(1<<fdmt->order >= fdmt->nf);
	assert(nbeams >= 1);

	// TODO: CHeck it's important that fmin < fmax??
	assert(fmin < fmax);
	assert(fmin > 0);
	assert(fmax > 0);

	//deltaT = int(np.ceil((maxDT-1) *(1./f_min**2 - 1./(f_min + deltaF)**2) / (1./f_min**2 - 1./f_max**2)))

	//fdmt->delta_t = (int)(ceilf((fdmt->maxDT-1) *(isquaref(fdmt->f_min) - isquaref(fdmt->f_min + fdmt->delta_f)) / (isquaref(f_min) - isquaref(f_max))));

	// Delta_t here is the number of time samples the maximum DM trajectory traverses
	// In the lowest channel. It is equivalent to the number of Dm trials you need to do
	// In the lowest channel to get out to the highest DM we asked for.
	fdmt->delta_t = calc_delta_t(fdmt, fdmt->fmin, fdmt->fmin + fdmt->df);
	fdmt->delta_t += 1; // Slightly different definition to origiinal

	// Allocate states as ping-pong buffer
	fdmt->state_size = fdmt->nbeams * fdmt->nf*fdmt->delta_t * fdmt->max_dt;
	fdmt->state_nbytes = fdmt->state_size * sizeof(fdmt_dtype);
	for (int s = 0; s < 2; s++) {
		fdmt->states[s].nw = fdmt->nbeams;
		fdmt->states[s].nx = fdmt->nf;
		fdmt->states[s].ny = fdmt->delta_t;
		fdmt->states[s].nz = fdmt->max_dt;
		array4d_malloc(&fdmt->states[s]);
	}

	fdmt->ostate.nw = fdmt->nbeams;
	fdmt->ostate.nx = 1;
	fdmt->ostate.ny = fdmt->max_dt;
	fdmt->ostate.nz = fdmt->max_dt;
	array4d_malloc(&fdmt->ostate);
	array4d_cuda_memset(&fdmt->ostate, 0);

	// save iteration setup
	int s = 0;
	for (int iiter = 1; iiter < fdmt->order+1; iiter++) {
		array4d_t* curr_state = &fdmt->states[s];
		s = (s + 1) % 2;
		array4d_t* new_state = &fdmt->states[s];
		fdmt_save_iteration(fdmt, iiter, curr_state, new_state);
	}

	fdmt->execute_count = 0;

	return 0;
}

int fdmt_initialise(const fdmt_t* fdmt, const array3d_t* indata, array4d_t* state)
{

	// indata is 3D array: (nbeams, nf, nt)
	// State is a 4D array: (nbeams, nf, deltat, nt) ( for the moment)

	assert(indata->nx == fdmt->nbeams);
	assert(indata->ny == fdmt->nf);
	assert(indata->nz == fdmt->nt);

	state->nw = fdmt->nbeams;
	state->nx = fdmt->nf;
	state->ny = fdmt->delta_t;
	state->nz = fdmt->max_dt;

	// zero off the state
	bzero(state->d, state->nw*state->nx*state->ny*state->nz*sizeof(fdmt_dtype));
	// Assign initial data to the state at delta_t=0
	for(int beam = 0 ; beam < fdmt->nbeams; beam++) {
		for (int c = 0; c < fdmt->nf; c++) {
			int outidx = array4d_idx(state, beam, c, 0, 0);
			int inidx = array3d_idx(indata, beam, c, 0);
			for (int t = 0; t < fdmt->nt; t++) {
				state->d[outidx + t] = indata->d[inidx + t];
			}
		}
	}

	// do partial sums initialisation (Equation 20.)
	// This (like everything barak does) is done as a recursive sum

	for(int beam = 0; beam < fdmt->nbeams; beam++) {
		// For each frequency channel
		for (int c = 0; c < fdmt->nf; c++) {
			// For each delta_t, i.e. each single-channel DM trial
			for (int idt = 1; idt < fdmt->delta_t; idt++) {
				int outidx = array4d_idx(state, beam, c, idt, 0);
				int iidx = array4d_idx(state, beam, c, idt-1, 0);
				int imidx = array3d_idx(indata, beam, c, indata->nz - 1);

				// The state for dt=d = the state for dt=(d-1) + the time-reversed input sample
				// for each time
				// (TODO: Not including a missing overlap here)
				for (int j = idt; j < fdmt->nt; j++) {
					state->d[outidx + j] = state->d[iidx + j] + indata->d[imidx - j];
				}
			}
		}
	}

	return 0;

}

int fdmt_iteration(const fdmt_t* fdmt,
		const int iteration_num,
		const array4d_t* indata,
		array4d_t* outdata)
{
	float df = fdmt->df; // channel resolution
	float delta_f = (float)(1 << iteration_num) * df; // Resolution of current iteration
	int delta_t = calc_delta_t(fdmt, fdmt->fmin, fdmt->fmin+delta_f); // Max DM

	// Outdata has size (nbeams, o_nf, o_nd1, fdmt->nt)
	outdata->nw = indata->nw;
	outdata->nx = indata->nx/2 + indata->nx % 2; // Add 1 to the frequency dimension if it's not divisible by 2
	outdata->ny = delta_t + 1;
	outdata->nz = indata->nz;

	assert(array4d_size(outdata) <= fdmt->state_size);

	//    printf("iteration %d df %f delta_f %f delta_t %d output nx=%d ny=%d nz%d\n",
	//           iteration_num, df, delta_f, delta_t, outdata->nx, outdata->ny, outdata->nz);

	// zero that output baby
	bzero(outdata->d, outdata->nw*outdata->nx * outdata->ny * outdata->nz * sizeof(fdmt_dtype));
	array4d_cuda_memset(outdata, 0);

	int shift_input = 0; // ?
	int shift_output = 0; // ?

	float fjumps = (float)outdata->nx; // Output number of channels
	float frange = fdmt->fmax - fdmt->fmin; // Width of band
	float fmin = fdmt->fmin; // Bottom of band

	float correction = 0.0;
	if (iteration_num > 0) {
		correction = df/2.0;
	}

	assert(indata->nw == fdmt->nbeams);
	// For each output sub-band
	for (int iif = 0; iif < outdata->nx; iif++) {
		float f_start = frange/fjumps * (float)iif + fmin; // Top freq of subband
		float f_end = frange/fjumps*((float)iif + 1) + fmin; // Bottom freq of subband
		float f_middle = (f_end - f_start)/2.0 + f_start - correction; // Middle freq of subband, less 0.5xresolution
		float f_middle_larger = (f_end - f_start)/2.0 + f_start + correction; // Middle freq of subband + 0.5x resolution (helps with rounding)

		// Max DM for this subband
		int delta_t_local = calc_delta_t(fdmt, f_start, f_end) + 1;

		// For each DM relevant for this subband
		for (int idt = 0; idt < delta_t_local; idt++) {
			int dt_middle = roundf(idt * cff(f_middle, f_start, f_end, f_start)); // Dt for middle freq less 0.5xresolution
			int dt_middle_index = dt_middle + shift_input;
			int dt_middle_larger = roundf(idt * cff(f_middle_larger, f_start, f_end, f_start)); // Dt for middle freq +0.5x resolution
			int dt_rest = idt - dt_middle_larger;
			int dt_rest_index = dt_rest + shift_input;

			int itmin = 0;
			int itmax = dt_middle_larger;

			//Output[i_F,i_dT + ShiftOutput,i_T_min:i_T_max] = Input[2*i_F, dT_middle_index,i_T_min:i_T_max];
			//int outidx = array4d_idx(outdata, beam, iif, idt+shift_output, 0);
			//int inidx1  = array4d_idx(indata, beam, 2*iif, dt_middle_index, 0);


			//printf("iteration %d channel %d freq %f idt %d dt_local "
			//	  "%d dt_middle %d dt_middle_larger %d dt_rest %d\n",
			//	  iteration_num, iif, f_middle, idt, delta_t_local, dt_middle_index, dt_middle_larger, dt_rest_index);

			// Here we handle the edge effects and set
			// OUtput state[freq, idx, 0:dtmin] = input_state[2xfreq, dt_middle, 0:dtmin]
			// where the DM would have overun the available times
			// This needs to be fixed for more careful time overlapping
			coord3_t dst_start = {.x = iif, .y = idt+shift_output, .z = 0};
			coord3_t src1_start = {.x = 2*iif, .y = dt_middle_index, .z = 0};

			// NB: Thisis run wit zcounts of 0, 1 and 2 - which break
			array_gpu_copy1(outdata, indata, &dst_start, &src1_start, dt_middle_larger);
			//cpu_copy2(&outdata->d[outidx + itmin], &indata->d[inidx1 + itmin], (itmax - itmin));
			//for (int i = itmin; i < itmax; i++) {
			//outdata->d[outidx + i] = indata->d[inidx1 + i];
			//}


			// Now we work on the remaining times that are guaranteed not to overrun the input dimensions
			itmin = dt_middle_larger;
			itmax = fdmt->max_dt;

			coord3_t src2_start = {.x = 2*iif + 1, .y = dt_rest_index, .z = 0};
			// src and dst now start from a bit offset
			src1_start.z = dt_middle_larger;
			dst_start.z = dt_middle_larger;
			int zcount = itmax - itmin;

			int maxt = dt_middle_larger;
			int src1_offset = array4d_idx(indata, 0, 2*iif, dt_middle_index, 0);
			int src2_offset = array4d_idx(indata, 0, 2*iif+1, dt_rest_index, 0);
			int out_offset = array4d_idx(outdata, 0, iif, idt, 0);
			//printf("iter %d iif %03d idt %02d src1_off %06d src2_off %06d out_off %06d maxt %02d dtmid %d dtr %d dtmidlg %d in [%d,%d,%d,%d]\n",
			//		iteration_num, iif, idt, src1_offset, src2_offset, out_offset, maxt, dt_middle_index, dt_rest_index, dt_middle_larger, indata->nw, indata->nx, indata->ny, indata->nz);


			if (2*iif + 1 < indata->nx) { // If the input data has this channel, we'll add it in
				//Output[i_F,i_dT + ShiftOutput,i_T_min:i_T_max] = Input[2*i_F, dT_middle_index,i_T_min:i_T_max] + Input[2*i_F+1, dT_rest_index,i_T_min - dT_middle_larger:i_T_max-dT_middle_larger]
				// playinga trick here - we're always addign the fastest moving index
				// Putting -dt_middle_larger in array3d_idx would have caused an assertion failure
				// But ofsetting by dt_middle_larger at the end, we get the best of all worlds

				//int inidx2 = array4d_idx(indata, beam, 2*iif+1, dt_rest_index, 0) - dt_middle_larger;

				array_gpu_sum1(outdata, indata, &dst_start, &src1_start, &src2_start, zcount);

				//for(int i = itmin; i < itmax; i++) {
				//  outdata->d[outidx + i] = indata->d[inidx1 + i] + indata->d[inidx2 + i];
				//}
				//cpu_sum1(&outdata->d[outidx + itmin], &indata->d[inidx1+itmin], &indata->d[inidx2+itmin], itmax-itmin);


			} else { // Just copy the input over. which basically assumes the upper channel is flaggedd/0
				// TODO: Could probably be done outside the iif loop to save evalutating IFs, but
				// Too tricky for the moment.
				//cpu_copy2(&outdata->d[outidx + itmin], &indata->d[inidx1 + itmin], (itmax - itmin));
				//for(int i = itmin; i < itmax; i++) {
				//  outdata->d[outidx + i] = indata->d[inidx1 + i];
				//	}
				array_gpu_copy1(outdata, indata, &dst_start, &src1_start, zcount);
			}
		}

	}
	return 0;
}


void __global__ cuda_fdmt_iteration_kernel2(float fmin, float frange, float fjumps, float correction)
{
	int beamno = blockIdx.x;
	int iif = blockIdx.y;
	int max_dt = blockDim.x;
	float f_start = frange/fjumps * (float)iif + fmin; // Top freq of subband
	float f_end = frange/fjumps*((float)iif + 1) + fmin; // Bottom freq of subband
	float f_middle = (f_end - f_start)/2.0 + f_start - correction; // Middle freq of subband, less 0.5xresolution
	float f_middle_larger = (f_end - f_start)/2.0 + f_start + correction; // Middle freq of subband + 0.5x resolution (helps with rounding)

	// Max DM for this subband
	//int delta_t_local = calc_delta_t(fdmt, f_start, f_end) + 1;
	int delta_t_local = 7; // TODO: Calculate
	int shift_input = 0;
	int shift_output = 0;


	// For each DM relevant for this subband
	for (int idt = 0; idt < delta_t_local; idt++) {
		int dt_middle = roundf(idt * cff(f_middle, f_start, f_end, f_start)); // Dt for middle freq less 0.5xresolution
		int dt_middle_index = dt_middle + shift_input;
		int dt_middle_larger = roundf(idt * cff(f_middle_larger, f_start, f_end, f_start)); // Dt for middle freq +0.5x resolution
		int dt_rest = idt - dt_middle_larger;
		int dt_rest_index = dt_rest + shift_input;

		int itmin = 0;
		int itmax = dt_middle_larger;

		//Output[i_F,i_dT + ShiftOutput,i_T_min:i_T_max] = Input[2*i_F, dT_middle_index,i_T_min:i_T_max];
		//int outidx = array4d_idx(outdata, beam, iif, idt+shift_output, 0);
		//int inidx1  = array4d_idx(indata, beam, 2*iif, dt_middle_index, 0);


		//printf("iteration %d channel %d freq %f idt %d dt_local "
		//	  "%d dt_middle %d dt_middle_larger %d dt_rest %d\n",
		//	  iteration_num, iif, f_middle, idt, delta_t_local, dt_middle_index, dt_middle_larger, dt_rest_index);

		// Here we handle the edge effects and set
		// OUtput state[freq, idx, 0:dtmin] = input_state[2xfreq, dt_middle, 0:dtmin]
		// where the DM would have overun the available times
		// This needs to be fixed for more careful time overlapping
		coord3_t dst_start = {.x = iif, .y = idt+shift_output, .z = 0};
		coord3_t src1_start = {.x = 2*iif, .y = dt_middle_index, .z = 0};
		//		array_gpu_copy1(outdata, indata, &dst_start, &src1_start, dt_middle_larger);
		//cpu_copy2(&outdata->d[outidx + itmin], &indata->d[inidx1 + itmin], (itmax - itmin));
		//for (int i = itmin; i < itmax; i++) {
		//outdata->d[outidx + i] = indata->d[inidx1 + i];
		//}


		// Now we work on the remaining times that are guaranteed not to overrun the input dimensions
		itmin = dt_middle_larger;
		itmax = max_dt;

		coord3_t src2_start = {.x = 2*iif + 1, .y = dt_rest_index, .z = 0};
		// src and dst now start from a bit offset
		src1_start.z = dt_middle_larger;
		dst_start.z = dt_middle_larger;
		int zcount = itmax - itmin;


		//		if (2*iif + 1 < indata->nx) { // If the input data has this channel, we'll add it in
		//			//Output[i_F,i_dT + ShiftOutput,i_T_min:i_T_max] = Input[2*i_F, dT_middle_index,i_T_min:i_T_max] + Input[2*i_F+1, dT_rest_index,i_T_min - dT_middle_larger:i_T_max-dT_middle_larger]
		//			// playinga trick here - we're always addign the fastest moving index
		//			// Putting -dt_middle_larger in array3d_idx would have caused an assertion failure
		//			// But ofsetting by dt_middle_larger at the end, we get the best of all worlds
		//
		//			//int inidx2 = array4d_idx(indata, beam, 2*iif+1, dt_rest_index, 0) - dt_middle_larger;
		//
		//			//array_gpu_sum1(outdata, indata, &dst_start, &src1_start, &src2_start, zcount);
		//
		//			//for(int i = itmin; i < itmax; i++) {
		//			//  outdata->d[outidx + i] = indata->d[inidx1 + i] + indata->d[inidx2 + i];
		//			//}
		//			//cpu_sum1(&outdata->d[outidx + itmin], &indata->d[inidx1+itmin], &indata->d[inidx2+itmin], itmax-itmin);
		//
		//
		//		} else { // Just copy the input over. which basically assumes the upper channel is flaggedd/0
		//			// TODO: Could probably be done outside the iif loop to save evalutating IFs, but
		//			// Too tricky for the moment.
		//			//cpu_copy2(&outdata->d[outidx + itmin], &indata->d[inidx1 + itmin], (itmax - itmin));
		//			//for(int i = itmin; i < itmax; i++) {
		//			//  outdata->d[outidx + i] = indata->d[inidx1 + i];
		//			//	}
		//			//array_gpu_copy1(outdata, indata, &dst_start, &src1_start, zcount);
		//		}
	}
}

// This version is totally hopeless - apparently it uses 44 register/thread, which makes it painfully slow, as the occupancy is low.
// launch bonds makes no difference. Gee it's so painfully crap, I'm kindof embarrassed.
void __global__
cuda_fdmt_iteration_kernel3(float fmin, float frange, float fjumps, float correction)
{
	int beamno = blockIdx.x;
	int iif = blockIdx.y;
	int max_dt = blockDim.x;
	float f_start = frange/fjumps * (float)iif + fmin; // Top freq of subband
	float f_end = frange/fjumps*((float)iif + 1) + fmin; // Bottom freq of subband
	float f_middle = (f_end - f_start)/2.0 + f_start - correction; // Middle freq of subband, less 0.5xresolution
	float f_middle_larger = (f_end - f_start)/2.0 + f_start + correction; // Middle freq of subband + 0.5x resolution (helps with rounding)
	//
	//	// Max DM for this subband
	//int delta_t_local = calc_delta_t(fdmt, f_start, f_end) + 1;
	int delta_t_local = 7; // TODO: Calculate
	int shift_input = 0;
	int shift_output = 0;
}

__host__ void cuda_fdmt_iteration3(const fdmt_t* fdmt, const int iteration_num, const array4d_t* indata, array4d_t* outdata)
{
	float df = fdmt->df; // channel resolution
	float delta_f = (float)(1 << iteration_num) * df; // Resolution of current iteration
	int delta_t = calc_delta_t(fdmt, fdmt->fmin, fdmt->fmin+delta_f); // Max DM

	// Outdata has size (nbeams, o_nf, o_nd1, fdmt->nt)
	outdata->nw = indata->nw;
	outdata->nx = indata->nx/2 + indata->nx % 2; // Add 1 to the frequency dimension if it's not divisible by 2
	outdata->ny = delta_t + 1;
	outdata->nz = indata->nz;

	assert(array4d_size(outdata) <= fdmt->state_size);


	//    printf("iteration %d df %f delta_f %f delta_t %d output nx=%d ny=%d nz%d\n",
	//           iteration_num, df, delta_f, delta_t, outdata->nx, outdata->ny, outdata->nz);

	// zero that output baby
	//bzero(outdata->d, outdata->nw*outdata->nx * outdata->ny * outdata->nz * sizeof(fdmt_dtype));
	array4d_cuda_memset(outdata, 0);
	float fjumps = (float)outdata->nx; // Output number of channels
	float frange = fdmt->fmax - fdmt->fmin; // Width of band
	float fmin = fdmt->fmin; // Bottom of band

	float correction = 0.0;
	if (iteration_num > 0) {
		correction = df/2.0;
	}

	assert(indata->nw == fdmt->nbeams);
	dim3 grid_shape(fdmt->nbeams, outdata->nx); // grid is (nbeams x nchannels wide)

	cuda_fdmt_iteration_kernel3<<<grid_shape, fdmt->max_dt>>>( fmin,  frange,  fjumps,  correction); // loops over all beams, subbands and dts
}

__global__ void cuda_fdmt_iteration_kernel4_sum (
		fdmt_dtype*  __restrict__ outdata,
		const fdmt_dtype* __restrict__ indata,
		int src_beam_stride,
		int dst_beam_stride,
		int delta_t_local,
		const int* __restrict__ ts_data)

{
	int beamno = blockIdx.x;
	int t = threadIdx.x;
	int nt = blockDim.x;
	fdmt_dtype* outp = outdata + beamno*dst_beam_stride + t;
	const fdmt_dtype* inp = indata + beamno*src_beam_stride + t;
	const int* ts_ptr = ts_data;
	if (t >= nt + delta_t_local) { // strictly could be the actuall delta_t for the out_offset, but who's counting?
		return;
	}
	for (int idt = 0; idt < delta_t_local; idt++ ) {
		int src1_offset = ts_ptr[0];
		int src2_offset = ts_ptr[1];
		int out_offset = ts_ptr[2];
		int mint = ts_ptr[3];
		if (t < mint) {
			outp[out_offset] = inp[src1_offset];
		} else {
			outp[out_offset] = inp[src1_offset] + inp[src2_offset];
		}

		ts_ptr += 4;
	}
}

__global__ void cuda_fdmt_iteration_kernel4_copy (
		fdmt_dtype*  __restrict__ outdata,
		const fdmt_dtype* __restrict__ indata,
		int src_beam_stride,
		int dst_beam_stride,
		int delta_t_local,
		const int* __restrict__ ts_data)

{
	int beamno = blockIdx.x;
	int t = threadIdx.x;
	fdmt_dtype* outp = outdata + beamno*dst_beam_stride + t;
	const fdmt_dtype* inp = indata + beamno*src_beam_stride + t;
	const int* ts_ptr = ts_data;
	for (int idt = 0; idt < delta_t_local; idt++ ) {
		int src1_offset = ts_ptr[0];
		int out_offset = ts_ptr[2];
		outp[out_offset] = inp[src1_offset];
		ts_ptr += 4;
	}
}

__host__ void cuda_fdmt_iteration4(const fdmt_t* fdmt, const int iteration_num, const array4d_t* indata, array4d_t* outdata)
{
	FdmtIteration* iter = fdmt->iterations.at(iteration_num-1);
	outdata->nw = iter->state_shape.w;
	outdata->nx = iter->state_shape.x;
	outdata->ny = iter->state_shape.y;
	outdata->nz = iter->state_shape.z;
	assert(indata->nz == outdata->nz);
	int nt = fdmt->nt;
	assert(array4d_size(outdata) <= fdmt->state_size);
	assert(outdata->nx == iter->dt_data.size());

	int src_beam_stride = array4d_idx(indata, 1, 0, 0, 0);
	int dst_beam_stride = array4d_idx(outdata, 1, 0, 0, 0);

	array4d_cuda_memset(outdata, 0);

	// For each output sub-band
	for (int iif = 0; iif < outdata->nx; iif++) {
		int* ts_data = iter->dt_data.at(iif)->d_device;
		int delta_t_local = iter->delta_ts.at(iif);
		const fdmt_dtype* src_start = &indata->d_device[0];
		fdmt_dtype* dst_start = &outdata->d_device[0];
		int tmax = min(nt + delta_t_local, indata->nz-1);
		//printf("iter %d iif %d TMAX %d nt  %d delta_t_local %d %d \n", iteration_num, iif, tmax, nt, delta_t_local, indata->nz);

		assert(tmax < indata->nz);


		// WARNING: This is sub-optimal as it doesn't use an integral number of warps, and
		// Can exceed teh maximum thread limit of the GPU , and use the threads sub-optimally.
		// More thought required to do this right.
		int nthreads = tmax;

		if(2*iif + 1 < indata->nx) { // do sum if there's a channel to sum
			cuda_fdmt_iteration_kernel4_sum<<<fdmt->nbeams, nthreads>>>(dst_start, src_start,
					src_beam_stride,
					dst_beam_stride,
					delta_t_local, ts_data);
			gpuErrchk(cudaPeekAtLastError());
		} else { // Do copy if there's no channel to add
			cuda_fdmt_iteration_kernel4_copy<<<fdmt->nbeams, nthreads>>>(dst_start, src_start,
					src_beam_stride,
					dst_beam_stride,
					delta_t_local, ts_data);
			gpuErrchk(cudaPeekAtLastError());
		}

	}


}

__global__ void cuda_fdmt_update_ostate(fdmt_dtype* __restrict__ ostate,
										const fdmt_dtype* __restrict__ indata)
{
	// Adds the indata into the ostate and shifts the ostate where the's some to add
	// shape of both arrays is [nbeams, max_dt, max_dt]
	// The time axis is the last one (it just has a big shape)
	// Where nbeams = gridDim.x
	// max_dt = gridDim.y
	// nt = blockDim.x = number of threads
	// Elsewhere we check that max_dt is a multiple of nt
	int ibeam = blockIdx.x;
	int nbeams = gridDim.x;

	int idt = blockIdx.y;
	int max_dt = gridDim.y;

	int nt = blockDim.x;
	int off = max_dt*(idt + ibeam*nbeams);
	int t = threadIdx.x;
	fdmt_dtype* optr = ostate + off;
	const fdmt_dtype* iptr = indata + off;

	// Add the new state for all but the last block
	while (t < max_dt - nt) {
		optr[t] = iptr[t] + optr[t + nt];

		// sync threads before doing the next block otherwise we don't copy the ostate correctly
		__syncthreads();
		t += nt;
	}
	// just copy the last block
	optr[t] = iptr[t];
}


__host__ void fdmt_update_ostate(fdmt_t* fdmt)
{
	// Run this after fdmt_execute_iterations when you want to take the output of the FDMT and update
	// the output state . i.e. do the delay and sum operation
	assert(fdmt->max_dt % fdmt->nt == 0);
	int s = fdmt->curr_state_idx;
	array4d_t* currstate = &fdmt->states[s];
	printf("w %d x %d y %d z %d\n", currstate->nw, currstate->nx, currstate->ny , currstate->nz);
	assert(currstate->nw == fdmt->ostate.nw);
	assert(currstate->nx == fdmt->ostate.nx);
	assert(currstate->ny == fdmt->ostate.ny);
	assert(currstate->nz == fdmt->ostate.nz);

	dim3 grid_shape(fdmt->nbeams, fdmt->max_dt);
	cuda_fdmt_update_ostate<<<grid_shape, fdmt->nt>>>(fdmt->ostate.d_device,
			currstate->d_device);

	gpuErrchk(cudaDeviceSynchronize());
}

__host__ void fdmt_copy_valid_ostate(fdmt_t* fdmt, fdmt_dtype* out)
{ // Copies only the first nt elemtns of the ostate to the output
	for(int b = 0; b < fdmt->nbeams; ++b) {
		for(int idt = 0; idt < fdmt->max_dt; ++idt) {
			int idx = array4d_idx(&fdmt->ostate, b, 0, idt, 0);
			gpuErrchk(cudaMemcpy(out + idx, fdmt->ostate.d_device + idx, sizeof(fdmt_dtype)*fdmt->nt, cudaMemcpyDeviceToHost));
		}
	}

}

__host__ void fdmt_copy_valid_ostate3(fdmt_t* fdmt, fdmt_dtype* out)
{ // Copies only the first nt elemtns of the ostate to the output
	bzero(out, sizeof(fdmt_dtype)*fdmt->nbeams*fdmt->max_dt*fdmt->nt);
	array4d_copy_to_host(&fdmt->ostate);
	int i = 0;
	for (int b = 0; b < fdmt->nbeams; ++b) {
		for(int idt = 0; idt < fdmt->max_dt; ++idt) {
			for (int t = 0; t < fdmt->nt; ++t) {
				int inidx = array4d_idx(&fdmt->ostate, b, 0, idt, t);
				int outidx = array4d_idx(&fdmt->ostate, b, 0, t, idt);
				out[outidx] = fdmt->ostate.d[inidx];
				i += 1;
			}
		}
	}
}

__host__ void fdmt_copy_valid_ostate2(const fdmt_t* fdmt, fdmt_dtype* dst)
{ // Copies only the first nt elemtns of the ostate to the output
	size_t dpitch = sizeof(fdmt_dtype)*fdmt->nt;
	size_t spitch = sizeof(fdmt_dtype)*fdmt->max_dt;
	size_t width = sizeof(fdmt_dtype)*fdmt->nt;
	size_t height = fdmt->max_dt;
	for(int b = 0; b < fdmt->nbeams; ++b) {
		int idx = array4d_idx(&fdmt->ostate, b, 0, 0, 0);
		gpuErrchk(cudaMemcpy2DAsync(&dst[idx], dpitch,
				&fdmt->ostate.d_device[idx], spitch,
				width, height, cudaMemcpyDeviceToHost));
	}
	gpuErrchk(cudaDeviceSynchronize());

}

int fdmt_execute_iterations(fdmt_t* fdmt)
{
	// Assumes data have been initialised into state[0]

	// Start that puppy up
	int s = 0;
	CudaTimer t;
	t.start();
	for (int iter = 1; iter < fdmt->order+1; iter++) {
		array4d_t* currstate = &fdmt->states[s];
		s = (s + 1) % 2;
		array4d_t* newstate = &fdmt->states[s];

		//fdmt_iteration(fdmt, iter, currstate, newstate);

		cuda_fdmt_iteration4(fdmt, iter, currstate, newstate);
		gpuErrchk(cudaPeekAtLastError());
		gpuErrchk(cudaDeviceSynchronize());
#ifdef DUMP_STATE
		char buf[128];
		array4d_copy_to_host(newstate);
		sprintf(buf, "state_s%d.dat", iter);
		array4d_dump(newstate, buf);
#endif
		//printf("Finisehd iteration %d\n", iter);

	}
	t.stop();
	cout << "FDMT Iterations only took " << t << endl;
	fdmt->curr_state_idx = s; // Tell people where to find the current state
}

int fdmt_execute(fdmt_t* fdmt, fdmt_dtype* indata, fdmt_dtype* outdata)
{
	array3d_t inarr = {.nx = fdmt->nbeams, .ny = fdmt->nf, .nz = fdmt->nt};
	inarr.d = indata;

	// Make the final outstate - this saves a memcpy on the final iteration
	array4d_t outstate;
	outstate.nw = fdmt->ostate.nw;
	outstate.nx = fdmt->ostate.nx;
	outstate.ny = fdmt->ostate.ny;
	outstate.nz = fdmt->ostate.nz;
	outstate.d = outdata;
	outstate.d_device = fdmt->ostate.d_device;

	// Initialise state
	CudaTimer tinit;
	tinit.start();
	int s = 0;
	fdmt_initialise(fdmt, &inarr, &fdmt->states[s]);
	tinit.stop();
	cout << "Initialisation took " << tinit << endl;

#ifdef DUMP_STATE
	// dump init state to disk
	char buf[128];
	sprintf(buf, "state_s%d.dat", 0);
	array4d_dump(&fdmt->states[s], buf);
	sprintf(buf, "initstate_e%d.dat", fdmt->execute_count);
	array4d_dump(&fdmt->states[s], buf);
#endif

	// copy initialised data to device
	CudaTimer tc;
	tc.start();
	array4d_copy_to_device(&fdmt->states[s]);
	tc.stop();
	cout << "Copy took " << tc << endl;

	// actually execute the iterations on the GPU
	fdmt_execute_iterations(fdmt);

	CudaTimer tupdate;
	tupdate.start();
	fdmt_update_ostate(fdmt);
	tupdate.stop();
	cout << "Delay and sum update took " << tupdate << endl;


#ifdef DUMP_STATE
	sprintf(buf, "ostate_e%d.dat", fdmt->execute_count);
	array4d_copy_to_host(&fdmt->ostate);
	array4d_dump(&fdmt->ostate, buf);
	array4d_t* currstate = &fdmt->states[fdmt->curr_state_idx];
	sprintf(buf, "finalstate_e%d.dat", fdmt->execute_count);
	array4d_copy_to_host(currstate);
	array4d_dump(currstate, buf);
#endif

	CudaTimer tback;
	tback.start();
	fdmt_copy_valid_ostate3(fdmt, outdata);
	tback.stop();
	cout << "Copy back to host took " << tback << endl;

	//printf("Returing form execute\n");

	fdmt->execute_count += 1;

	return 0;
}
