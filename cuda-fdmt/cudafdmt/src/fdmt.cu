#include <iostream>
#include "cpu_kernels.h"
#include "gpu_kernels.h"
#include "fdmt.h"
#include "fdmt_utils.h"
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

__host__ __device__ int calc_delta_t(float f_start, float f_end, float fmin, float fmax, int max_dt)
{
	float rf = cff(f_start, f_end, fmin, fmax);
	float delta_tf = ((float)max_dt-1.0) * rf;
	int delta_t = (int)ceilf(delta_tf);

	return delta_t;
}

__host__ FdmtIteration* fdmt_save_iteration(fdmt_t* fdmt, const int iteration_num, const array4d_t* indata, array4d_t* outdata)
{

	// If the number of incoming channels is not a multiple of two, we copy the highest channel to the next state
	// without doing anything to it (except making sure it lives in the right place in the destination state.
	// We track the frequency resolution of the top channel separately from all the other bottom channels

	const float df = fdmt->df; // channel resolution
	const float fmin = fdmt->fmin; // Bottom of band

	const bool do_copy = indata->nx %2 == 1; // true if the top channel will be copied rather than FDMTed
	const int nf = indata->nx/2 + indata->nx % 2; 	// Number of channels in the output state. Adds 1 if the top channel will be copied.

	if (do_copy) {
		fdmt->_df_top += 0; //top channel width unchanged
	} else {
		fdmt->_df_top += fdmt->_df_bot; // We will be wider by a new channel
	}

	fdmt->_df_bot *= 2.0; // Bottom channels will be added together

	float delta_f;
	if (nf == 1) { // This is the last iteration
		delta_f = fdmt->_df_top;
	} else {
		delta_f = fdmt->_df_bot;
	}

	const float fres = fdmt->_df_bot;
	int delta_t = calc_delta_t(fdmt, fdmt->fmin, fdmt->fmin+delta_f); // Max DM
	int ndt = delta_t + 1;

	// Outdata has size (nbeams, o_nf, o_nd1, fdmt->nt)
	outdata->nw = indata->nw;
	outdata->nx = nf;
	outdata->ny = ndt;
	outdata->nz = fdmt->nt + ndt; // everything else is zeros

	//assert(array4d_size(outdata) <= fdmt->state_size);

	FdmtIteration* iter = new FdmtIteration(outdata->nw, outdata->nx, outdata->ny, outdata->nz);
	fdmt->iterations.push_back(iter);


//	printf("Iteration %d max_dt %d nf=%d fres=%f fmin=%f inshape: [%d, %d, %d, %d] outshape: [%d, %d, %d, %d]\n",
//			iteration_num, fdmt->max_dt, nf, fres, fmin,
//				indata->nw, indata->nx, indata->ny, indata->nz,
//				outdata->nw, outdata->nx, outdata->ny, outdata->nz);

	float correction = 0.0;
	if (iteration_num > 0) {
		correction = df/2.0;
	}

	int shift_input = 0;
	int shift_output = 0;

	// For each output sub-band
	for (int iif = 0; iif < nf; iif++) {
		bool copy_subband = false;
		float f_start = fres * (float)iif + fmin - df/2.0f; // Freq of bottom output subband
		float f_end;
		float f_middle;
		int delta_t_local;

		if (iif < nf - 1) { // all the bottom subbands
			f_end = f_start + fres;
			f_middle = f_start +  fres/2.0f - correction; // Middle freq of subband, less 0.5xresolution
			delta_t_local = calc_delta_t(fdmt, f_start, f_end) + 1;
		} else  { // if this is the top output subband
			assert(iif == nf - 1);
			//printf("iif %d final subband\n", iif);

			if (do_copy) { // this iteration is a copy iteration -  there is no subband above this one in the input data
				// the output channel width equals the input channel width
				//printf("no Subband available. Copy: fdmt->_df_top=%f\n", fdmt->_df_top);
				f_end = f_start + fdmt->_df_top*2.0;
				f_middle = f_start + fdmt->_df_top - correction; // Middle freq of subband, less 0.5xresolution
				// Tell that code down there to mark this subband to copy the output across
				copy_subband = true;
				delta_t_local = fdmt->_ndt_top;

			} else { // There are 2 subbands available in the input data
				// The width of the output subband is the sum of the input suband (which is fres/2.0)
				// plus whatever the previous output was
				f_end = f_start + fdmt->_df_top;
				f_middle = f_start + fres/2.0 - correction;
				//printf("2 subbands available: fdmt->_df_top=%f. fres %f f_start %f f_end %f\n", fdmt->_df_top, fres, f_start, f_end);
				delta_t_local = calc_delta_t(fdmt, f_start, f_end) + 1;
				fdmt->_ndt_top = delta_t_local;
			}
		}

		float f_middle_larger = f_middle + 2*correction; // Middle freq of subband + 0.5x resolution (helps with rounding)


		// Note; we must not overwrite the max ndt - doign too many down low will give us too much resolution
		// Up high.
		if (delta_t_local > ndt) {
			//printf("YUK! delta_t_local %d > ndt %d\n", delta_t_local, ndt);
			delta_t_local = ndt;
		}
		iter->add_subband(delta_t_local);
//		printf("oif %d iif1=%d iif2=%d dt_loc=%d f_start %f f_end %f f_middle %f f_middle_larger %f\n", iif,
//					2*iif, 2*iif+1, delta_t_local, f_start, f_end, f_middle, f_middle_larger);
		if (iif == 0) {
			assert(delta_t_local == ndt);// Should populate all delta_t in the lowest band!!!
		}

		assert(f_start < f_middle);
		assert(f_middle < f_middle_larger);
		assert(f_middle_larger < f_end);


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
			//printf("idt %d dtmid %d dt_mid_l %d dt_rest %d\n", idt, dt_middle_index, dt_middle_larger, dt_rest);

			dim3 dst_start(iif, idt+shift_output,0);
			dim3 src1_start(2*iif, dt_middle_index, 0);
			dim3 src2_start(2*iif + 1, dt_rest_index, 0);

			//array_gpu_copy1(outdata, indata, &dst_start, &src1_start, dt_middle_larger);

     		// Now we work on the remaining times that are guaranteed not to overrun the input dimensions
			itmin = dt_middle_larger;
			itmax = fdmt->max_dt;
			assert(itmax > itmin);

# 			// save number of operations the FDMT takes
                        fdmt->nops += (itmax - itmin)*fdmt->nbeams;
			// src and dst now start from a bit offset
			src1_start.z = dt_middle_larger;
			dst_start.z = dt_middle_larger;

			// We shouldn't overrun the incoming array
			if (dt_middle_index >= indata->ny) {
				printf("Bleagh! dt_middle_index exceeded input array. iteration_num %d oif=%d idt %d dt_middle_index %d ny %d delta_t_local %d\n", iteration_num, iif, idt, dt_middle_index, indata->ny, delta_t_local);
				dt_middle_index = indata->ny -1 ; // NOT SURE I SHOULD BE DOING THIS!
			}
			if (dt_rest_index >= indata->ny) {
				printf("Bleagh! dt_rest_index exceeded input array. iteration_num %d oif=%d idt %d dt_middle_index %d ny %d delta_t_local %d\n", iteration_num, iif, idt, dt_middle_index, indata->ny, delta_t_local);
				dt_rest_index = indata->ny - 1;
			}
			assert(dt_middle_index < indata->ny);
			assert(dt_rest_index < indata->ny);
			// TODO: ADD MORE BOUNDS CHECKS

			int mint = dt_middle_larger;

			int src1_offset = array4d_idx(indata, 0, 2*iif, dt_middle_index, 0);
			int src2_offset = array4d_idx(indata, 0, 2*iif+1, dt_rest_index, 0) - mint; // src2 offset is by mint because it only starts adding after mint samples
			int out_offset = array4d_idx(outdata, 0, iif, idt, 0);
			if (copy_subband) {
				src1_offset = array4d_idx(indata, 0, 2*iif, idt, 0);
				src2_offset = -1;
				mint = 0;
			}
			//printf("Subband idt %d src1_offset %d src2_offset %d out_offset %d mint %d\n", idt, src1_offset, src2_offset, out_offset, mint);
			iter->save_subband_values(idt, src1_offset, src2_offset, out_offset, mint);
		}
	}
	iter->copy_to_device();

	return iter;
}



int fdmt_create(fdmt_t* fdmt, float fmin, float fmax, int nf, int max_dt, int nt, int nbeams, int nbeams_alloc, bool dump_data)
{
	// Expect center frequencies on the input here. Interanlly use the bottom edge frequency
	// nbeams_alloc = number of beams to allocate working memory for. <= nbeams. Beams will be processed serially
	// if <= 0, then process all beams
	fdmt->max_dt = max_dt;
	fdmt->nt = nt;
	fdmt->nf = nf;
	fdmt->df = (fmax - fmin)/((float) fdmt->nf);
	fdmt->fmin = fmin;
	fdmt->fmax = fmax;
	fdmt->order = (int)ceil(log(fdmt->nf)/log(2.0));
	fdmt->nbeams = nbeams;
	if (nbeams_alloc <= 0) {
		fdmt->nbeams_alloc = nbeams;
	} else {
		fdmt->nbeams_alloc = nbeams_alloc;
	}
	fdmt->dump_data = dump_data;
	fdmt->nops = 0;
	bool host_alloc = dump_data;
	assert(nf > 0);
	assert(max_dt > 0);
	assert(nt > 0);
	assert(fdmt->max_dt >= fdmt->nt);
	assert(fdmt->max_dt % fdmt->nt == 0); // max_dt needs to be a multipel of nt
	assert(1<<fdmt->order >= fdmt->nf);
	assert(nbeams >= 1);
	assert(fdmt->nbeams_alloc <= nbeams);
	assert(fdmt->nbeams_alloc > 0);

	// TODO: CHeck it's important that fmin < fmax??
	assert(fmin < fmax);
	assert(fmin > 0);
	assert(fmax > 0);

	//deltaT = int(np.ceil((maxDT-1) *(1./f_min**2 - 1./(f_min + deltaF)**2) / (1./f_min**2 - 1./f_max**2)))

	//fdmt->delta_t = (int)(ceilf((fdmt->maxDT-1) *(isquaref(fdmt->f_min) - isquaref(fdmt->f_min + fdmt->delta_f)) / (isquaref(f_min) - isquaref(f_max))));

	// Delta_t here is the number of time samples the maximum DM trajectory traverses
	// In the lowest channel. It is equivalent to the number of Dm trials you need to do
	// In the lowest channel to get out to the highest DM we asked for.
	fdmt->delta_t = calc_delta_t(fdmt, fdmt->fmin-fdmt->df/2., fdmt->fmin + fdmt->df/2.);
	fdmt->delta_t += 1; // Slightly different definition to original

	// Allocate states as ping-pong buffer
	for (int s = 0; s < 2; s++) {
		fdmt->states[s].nw = fdmt->nbeams_alloc;
		fdmt->states[s].nx = fdmt->nf;
		fdmt->states[s].ny = fdmt->delta_t;
		fdmt->states[s].nz = fdmt->delta_t + nt;
		//array4d_malloc(&fdmt->states[s], host_alloc, true);
	}

	// save iteration setup
	int s = 0;
	fdmt->_df_bot = fdmt->df;
	fdmt->_df_top = fdmt->df;
	array4d_t biggest_state = fdmt->states[s];
	size_t biggest_state_size = array4d_size(&biggest_state);

	for (int iiter = 1; iiter < fdmt->order+1; iiter++) {
		array4d_t* curr_state = &fdmt->states[s];
		s = (s + 1) % 2;
		array4d_t* new_state = &fdmt->states[s];
		fdmt_save_iteration(fdmt, iiter, curr_state, new_state);
		size_t new_state_size = array4d_size(new_state);
		if (new_state_size > biggest_state_size) {
			biggest_state = *new_state;
			biggest_state_size = new_state_size;
		}
	}

	printf("Biggest state is %d elements = %d MiB = ", biggest_state_size, biggest_state_size*sizeof(fdmt_dtype)/1024/1024);
	array4d_print_shape(&biggest_state);
	printf("\n");

	for(int s = 0; s < 2; ++s) {
		fdmt->states[s].nw = biggest_state.nw;
		fdmt->states[s].nx = biggest_state.nx;
		fdmt->states[s].ny = biggest_state.ny;
		fdmt->states[s].nz = biggest_state.nz;
		array4d_malloc(&fdmt->states[s], host_alloc, true);
	}

	fdmt->state_size = array4d_size(&fdmt->states[0]);
	fdmt->state_nbytes = fdmt->state_size * sizeof(fdmt_dtype);

	fdmt->ostate.nw = fdmt->nbeams;
	fdmt->ostate.nx = 1;
	fdmt->ostate.ny = fdmt->max_dt;
	fdmt->ostate.nz = fdmt->max_dt + nt;
	array4d_malloc(&fdmt->ostate, host_alloc, true);
	array4d_cuda_memset(&fdmt->ostate, 0);

	fdmt->weights.nw = 1;
	fdmt->weights.nx = 1;
	fdmt->weights.ny = 1;
	fdmt->weights.nz = fdmt->max_dt;
	array4d_malloc(&fdmt->weights, host_alloc, true);
	array4d_fill_device(&fdmt->weights, 1.0);



	fdmt->execute_count = 0;

	fdmt_calculate_weights(fdmt);

	return 0;
}

int fdmt_initialise(const fdmt_t* fdmt, const array3d_t* indata, array4d_t* state)
{

	// indata is 3D array: (nbeams, nf, nt)
	// State is a 4D array: (nbeams, nf, deltat, max_dt) ( for the moment)

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
					fdmt_dtype n = (fdmt_dtype)  j;
					state->d[outidx + j] = (state->d[iidx + j] + indata->d[imidx - j]);
					printf("Chan c=%c idt=%idt j=%d stateout=%f sattein=%f indata=%f\n", c, idt, j,
							state->d[outidx + j], state->d[iidx + j], indata->d[imidx - j]);
				}
			}
		}
	}

	return 0;

}

void __global__ fdmt_initialise_kernel(const fdmt_dtype* __restrict__ indata,
		fdmt_dtype* __restrict__ state, int delta_t, int max_dt, int nt)
{
	// indata is 4D array: (nbeams, nf, 1, nt): index [ibeam, c, 0, t] = t + nt*(0 + 1*(c + nf*ibeam))
	// State is a 4D array: (nbeams, nf, delta_t, max_dt) ( for the moment)
	// full index [ibeam, c, idt, t] is t + max_dt*(idt + delta_t*(c + nf*ibeam))

	int nbeams = gridDim.x; // number of beams
	int nf = blockDim.x; // Number of frequencies
	int ibeam = blockIdx.x; // beam number
	int c = threadIdx.x; // Channel number

	// Assign initial data to the state at delta_t=0
	int outidx = array4d_idx(nbeams, nf, delta_t, max_dt, ibeam, c, 0, 0);
	int imidx = array4d_idx(nbeams, nf, 1, nt, ibeam, c, 0, 0);
	for (int t = 0; t < nt; ++t) {
		state[outidx + t] = indata[imidx + t];
	}

	// Do partial sums initialisation recursively (Equation 20.)
	for (int idt = 1; idt < delta_t; ++idt) {
		int outidx = array4d_idx(nbeams, nf, delta_t, max_dt, ibeam, c, idt, 0);
		int iidx   = array4d_idx(nbeams, nf, delta_t, max_dt, ibeam, c, idt-1, 0);
		int imidx  = array4d_idx(nbeams, nf, 1, nt, ibeam, c, 0, nt -1 );

		// The state for dt=d = the state for dt=(d-1) + the time-reversed input sample
		// for each time
		// (TODO: Not including a missing overlap with the previous block here)
		// originally this was j=idt, rather than j=0. But that just meant that 0<=j<idt were zero, which seems weird.

		for(int j = 0; j < nt; ++j) {
			state[outidx + j] = (state[iidx + j] + indata[imidx -j]);
		}
	}
}

void __global__ fdmt_initialise_kernel2(const fdmt_dtype* __restrict__ indata,
		fdmt_dtype* __restrict__ state, int delta_t, int max_dt, int nt, bool count)
{
	// indata is 4D array: (nbeams, nf, 1, nt): index [ibeam, c, 0, t] = t + nt*(0 + 1*(c + nf*ibeam))
	// State is a 4D array: (nbeams, nf, delta_t, delta_t + nt) ( for the moment)
	// full index [ibeam, c, idt, t] is t + max_dt*(idt + delta_t*(c + nf*ibeam))
	// If coutn is true, it initialises the input to the number of cells that will be added

	int nbeams = gridDim.x; // number of beams
	int nf = gridDim.y; // Number of frequencies
	int tblock = blockDim.x; // number of samples per thread block
	int ibeam = blockIdx.x; // beam number
	int c = blockIdx.y; // Channel number
	int t = threadIdx.x; // sample number

	// Assign initial data to the state at delta_t=0
	int outidx = array4d_idx(nbeams, nf, delta_t, delta_t + nt, ibeam, c, 0, 0);
	int imidx = array4d_idx(nbeams, nf, 1, nt, ibeam, c, 0, 0);
	while (t < nt) {
		if (count) {
			state[outidx + t] = 1.;
		} else {
			state[outidx + t] = indata[imidx + t];
		}
		t += tblock;
	}

	// Do partial sums initialisation recursively (Equation 20.)
	for (int idt = 1; idt < delta_t; ++idt) {
		int outidx = array4d_idx(nbeams, nf, delta_t, delta_t + nt, ibeam, c, idt, 0);
		int iidx   = array4d_idx(nbeams, nf, delta_t, delta_t + nt, ibeam, c, idt-1, 0);
		int imidx  = array4d_idx(nbeams, nf, 1, nt, ibeam, c, 0, 0 );

		// The state for dt=d = the state for dt=(d-1) + the time-reversed input sample
		// for each time
		// (TODO: Not including a missing overlap with the previous block here)
		// originally this was j=idt, rather than j=0. But that just meant that 0<=j<idt were zero, which seems weird.
		t = threadIdx.x; // reset t
		fdmt_dtype c1 = (fdmt_dtype(idt));
		fdmt_dtype c2 = (fdmt_dtype(idt+1));
		//c1 = 1.;
		//c2 = 1.;
		while (t < nt) {
			if (count) {
				state[outidx + t] = fdmt_dtype(idt + 1);
			} else {
				state[outidx + t] = (state[iidx + t]*c1 + indata[imidx + t])/c2;
			}
			t += tblock;
		}
	}
}

int fdmt_initialise_gpu(const fdmt_t* fdmt, const array4d_t* indata, array4d_t* state, bool count)
{
	// indata is 4D array: (nbeams, nf, 1, nt)
	// State is a 4D array: (nbeams, nf, deltat, max_dt) ( for the moment)

	int nbeams = indata->nw; // number of beams in this batch
	assert(nbeams > 0);
	assert(nbeams <= fdmt->nbeams_alloc);

	if (! count) {
		assert(indata->nx == fdmt->nf);
		assert(indata->ny == 1);
		assert(indata->nz == fdmt->nt);
	}

	//state->nw = fdmt->nbeams;
	state->nw = nbeams; // nbeams in this batch
	state->nx = fdmt->nf;
	state->ny = fdmt->delta_t;
	state->nz = fdmt->nt + fdmt->delta_t;

	// zero off the state
	array4d_cuda_memset(state, 0);
	//gpuErrchk(cudaDeviceSynchronize());

	dim3 grid_shape(nbeams, fdmt->nf);
	//fdmt_initialise_kernel<<<fdmt->nbeams, fdmt->nf>>>(indata->d_device, state->d_device, fdmt->delta_t, fdmt->max_dt, fdmt->nt);
	int nthreads = 256;
	fdmt_initialise_kernel2<<<grid_shape, nthreads>>>(indata->d_device, state->d_device, fdmt->delta_t, fdmt->max_dt, fdmt->nt, count);
	//gpuErrchk(cudaDeviceSynchronize());

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


// The thing I like about this version is there' some small hope that it will cache some of teh data, because it loops over idt
// On the other hand, we probably should to it explicitly with shared memory if we want to do it prperly.
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

// The thing I like about this i that it has loads of blocks, so works well even if nbeams is small.
// But: it doesn't do much caching.
__global__ void cuda_fdmt_iteration_kernel5_sum (
		fdmt_dtype*  __restrict__ outdata,
		const fdmt_dtype* __restrict__ indata,
		int src_beam_stride,
		int dst_beam_stride,
		int tmax,
		int tend,
		const int* __restrict__ ts_data)

{
	int beamno = blockIdx.x;
	int idt = blockIdx.y;
	int ndt = gridDim.y;
	int t = threadIdx.x;
	int nt = blockDim.x;
	const int* ts_ptr = ts_data + 4*idt;

	int src1_offset = ts_ptr[0];
	int src2_offset = ts_ptr[1];
	int out_offset = ts_ptr[2];
	int mint = ts_ptr[3];

	fdmt_dtype* outp = outdata + beamno*dst_beam_stride + t;
	const fdmt_dtype* inp = indata + beamno*src_beam_stride + t;

	while(t < mint) {
		outp[out_offset] = inp[src1_offset];
		t += nt;
		outp += nt;
		inp += nt;
	}

	while(t < tmax) {
		outp[out_offset] = inp[src1_offset] + inp[src2_offset];
		t += nt;
		outp += nt;
		inp += nt;
	}

	int tend1 = min(tend, tmax + mint);

	while(t < tend1) {
		outp[out_offset] = inp[src2_offset];
		t += nt;
		outp += nt;
		inp += nt;
	}

	while(t < tend) {
		outp[out_offset] = 0;
		t += nt;
		outp += nt;
		inp += nt;
	}
}

// The thing I like about this i that it has loads of blocks, so works well even if nbeams is small.
// But: it doesn't do much caching.
__global__ void cuda_fdmt_iteration_kernel5_copy (
		fdmt_dtype*  __restrict__ outdata,
		const fdmt_dtype* __restrict__ indata,
		int src_beam_stride,
		int dst_beam_stride,
		int tend,
		const int* __restrict__ ts_data)

{
	int beamno = blockIdx.x;
	int idt = blockIdx.y;
	int t = threadIdx.x;
	int nt = blockDim.x;
	fdmt_dtype* outp = outdata + beamno*dst_beam_stride + t;
	const fdmt_dtype* inp = indata + beamno*src_beam_stride + t;
	const int* ts_ptr = ts_data + 4*idt;
	int src1_offset = ts_ptr[0];
	int out_offset = ts_ptr[2];

	while(t < tend) {
		outp[out_offset] = inp[src1_offset];
		t += nt;
		outp += nt;
		inp += nt;
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

__host__ void cuda_fdmt_iteration4(const fdmt_t* fdmt, const int iteration_num, const array4d_t* indata, array4d_t* outdata, int nbeams)
{
	FdmtIteration* iter = fdmt->iterations.at(iteration_num-1);
	//outdata->nw = iter->state_shape.w;
	outdata->nw = nbeams;
	outdata->nx = iter->state_shape.x;
	outdata->ny = iter->state_shape.y;
	outdata->nz = iter->state_shape.z;
	//assert(indata->nz == outdata->nz); // This almost certainly isn't the case - it grows with iteration

	int nt = fdmt->nt;
	assert(array4d_size(outdata) <= fdmt->state_size);
	assert(outdata->nx == iter->dt_data.size());

	int src_beam_stride = array4d_idx(indata, 1, 0, 0, 0);
	int dst_beam_stride = array4d_idx(outdata, 1, 0, 0, 0);

	// Not sure this memset is necessary.
	//array4d_cuda_memset(outdata, 0);

	// For each output sub-band
	for (int iif = 0; iif < outdata->nx; iif++) {
		int* ts_data = iter->dt_data.at(iif)->d_device;
		int delta_t_local = iter->delta_ts.at(iif);
		const fdmt_dtype* src_start = &indata->d_device[0];
		fdmt_dtype* dst_start = &outdata->d_device[0];
		int tmax = min(nt + delta_t_local, indata->nz);
		int tend = outdata->nz;
		assert(tmax <= indata->nz);
		assert(tmax <= tend);


		// WARNING: This is sub-optimal as it doesn't use an integral number of warps, and
		// Can exceed teh maximum thread limit of the GPU , and use the threads sub-optimally.
		// More thought required to do this right.
		// kernel4 requires nthreads = tax
		int nthreads = tmax;
		dim3 grid_size(nbeams, delta_t_local);

		//printf("Iteration %d iif %d indata->nz %d outdata->nz %d nt=%d delta_t_local %d tmax %d tend %d\n", iteration_num, iif,
				//indata->nz, outdata->nz, nt, delta_t_local, tmax, tend);


		int nthread = 128;
		if(2*iif + 1 < indata->nx) { // do sum if there's a channel to sum

			cuda_fdmt_iteration_kernel5_sum<<<grid_size, nthread>>>(dst_start, src_start,
					src_beam_stride,
					dst_beam_stride,
					tmax,tend,
					ts_data);
			//gpuErrchk(cudaPeekAtLastError());
		} else { // Do copy if there's no channel to add


			cuda_fdmt_iteration_kernel5_copy<<<grid_size, nthread>>>(dst_start, src_start,
								src_beam_stride,
								dst_beam_stride,
								tmax, ts_data);

			//gpuErrchk(cudaPeekAtLastError());
		}

	}


}

__global__ void cuda_fdmt_update_ostate(fdmt_dtype* __restrict__ ostate,
										const fdmt_dtype* __restrict__ indata,
										const fdmt_dtype weight,
										int nt)
{
	// Adds the indata into the ostate and shifts the ostate where the's some to add
	// shape of both arrays is [nbeams, max_dt, max_dt + nt]
	// The time axis is the last one (it just has a big shape)
	// Where nbeams = gridDim.x
	// max_dt = gridDim.y
	// nt = blockDim.x = number of threads
	// Elsewhere we check that max_dt is a multiple of nt

	int ibeam = blockIdx.x;
	int nbeams = gridDim.x;
	int blockt = blockDim.x;
	int max_dt = gridDim.y;
	int idt = blockIdx.y;

	int off = array4d_idx(1, nbeams, max_dt, max_dt+nt, 0, ibeam, idt, 0);
	fdmt_dtype* optr = ostate + off;
	const fdmt_dtype* iptr = indata + off;
	// Add the new state for all but the last block
	for (int t = threadIdx.x; t < max_dt + nt; t += blockt) {
		// makign this optr[t] = iptr[t+-1] + optr[t + nt] makes the DC RMS worse
		// optr[t] = iptr[t] + optr[t + nt +- 1]; also worse
		// So 		optr[t] = iptr[t] + optr[t + nt];
		// It's just weird that we have such a noisy response to DC input

		if (t < nt) {
			// Weight only the last block by the weights
			optr[t] = (iptr[t] + optr[t + nt])*weight;
		} else if (t >= max_dt) {
			optr[t] = iptr[t];
		} else {
			optr[t] = (iptr[t] + optr[t + nt]);
		}

		// sync threads before doing the next block otherwise we don't copy the ostate correctly
		__syncthreads();
	}
}


__host__ void fdmt_update_ostate(fdmt_t* fdmt, int ibeam, int nbeams)
{
	// Run this after fdmt_execute_iterations when you want to take the output of the FDMT and update
	// the output state . i.e. do the delay and sum operation
	assert(fdmt->max_dt % fdmt->nt == 0);
	int s = fdmt->curr_state_idx;
	array4d_t* currstate = &fdmt->states[s];
	//assert(currstate->nw == fdmt->ostate.nw);
	assert(currstate->nx == fdmt->ostate.nx);
	assert(currstate->ny == fdmt->ostate.ny);
	assert(currstate->nz == fdmt->ostate.nz);
	fdmt_dtype* optr = fdmt->ostate.d_device + array4d_idx(&fdmt->ostate, ibeam, 0, 0, 0);

	dim3 grid_shape(nbeams, fdmt->max_dt);
	cuda_fdmt_update_ostate<<<grid_shape, 256>>>(optr,
			currstate->d_device, rsqrtf(fdmt->nf), fdmt->nt);

	//gpuErrchk(cudaDeviceSynchronize());
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

__host__ void fdmt_copy_valid_ostate3(fdmt_t* fdmt, array4d_t* out)
{   // Copies only the first nt elemtns of the ostate to the output
	// Transposes data to beam, idt, t order
	// NOTE: THis is extremely inefficient! (it copies waay to much data over)
	// But: It works.
	array4d_copy_to_host(&fdmt->ostate);
	int i = 0;
	assert(out->nw == 1);
	assert(out->nx == fdmt->nbeams);
	assert(out->ny == fdmt->max_dt);
	assert(out->nz == fdmt->nt);

	for (int b = 0; b < fdmt->nbeams; ++b) {
		for(int idt = 0; idt < fdmt->max_dt; ++idt) {
			for (int t = 0; t < fdmt->nt; ++t) {
				int inidx = array4d_idx(&fdmt->ostate, b, 0, idt, t);
				int outidx = array4d_idx(out, 0, b, idt, t);
				out->d[outidx] = fdmt->ostate.d[inidx];
				i += 1;
			}
		}
	}
}

__host__ void fdmt_copy_valid_ostate2(const fdmt_t* fdmt, array4d_t* out)
{
	// ostate has size [nbeams, 1, max_dt, max_dt + nt]
	// dst has size [1, nbeams, max_dt, nt]
	assert(out->nw == 1);
	assert(out->nx == fdmt->nbeams);
	assert(out->ny == fdmt->max_dt);
	assert(out->nz == fdmt->nt);

	size_t dpitch = sizeof(fdmt_dtype)*fdmt->nt;
	size_t spitch = sizeof(fdmt_dtype)*(fdmt->max_dt + fdmt->nt);
	size_t width = sizeof(fdmt_dtype)*fdmt->nt;
	size_t height = fdmt->max_dt;
	for(int b = 0; b < fdmt->nbeams; ++b) {
		int inidx = array4d_idx(&fdmt->ostate, b, 0, 0, 0);
		int oidx = array4d_idx(out, 0, b, 0, 0);
		gpuErrchk(cudaMemcpy2D(&out->d[oidx], dpitch,
				&fdmt->ostate.d_device[inidx], spitch,
				width, height, cudaMemcpyDeviceToHost));
	}
}

int fdmt_execute_iterations(fdmt_t* fdmt, int nbeams)
{
	// Assumes data have been initialised into state[0]

	// Start that puppy up
	int s = 0;
	fdmt->t_iterations.start();
	for (int iter = 1; iter < fdmt->order+1; iter++) {
		array4d_t* currstate = &fdmt->states[s];
		s = (s + 1) % 2;
		array4d_t* newstate = &fdmt->states[s];
		cuda_fdmt_iteration4(fdmt, iter, currstate, newstate, nbeams);
		//gpuErrchk(cudaPeekAtLastError());
		//gpuErrchk(cudaDeviceSynchronize());
#ifdef DUMP_STATE
		char buf[128];
		array4d_copy_to_host(newstate);
		sprintf(buf, "state_s%d.dat", iter);
		array4d_dump(newstate, buf);
#endif
	}
	fdmt->t_iterations.stop();
	//cout << "FDMT Iterations only took " << t << endl;
	fdmt->curr_state_idx = s; // Tell people where to find the current state
}

int fdmt_execute_batch(fdmt_t* fdmt, fdmt_dtype* indata, fdmt_dtype* outdata, int ibeam, int nbeams)
{
	// Runs nbeams beams starting at ibeam through the FDMT
	// and updates the relevant output state
	assert(ibeam + nbeams <= fdmt->nbeams);

	array4d_t inarr;
	inarr.nw = nbeams;
	inarr.nx = fdmt->nf;
	inarr.ny = 1;
	inarr.nz = fdmt->nt;
	inarr.d_device = indata + ibeam*(fdmt->nt*fdmt->nf);

	// Initialise state
	fdmt->t_init.start();
	int s = 0;
	fdmt_initialise_gpu(fdmt, &inarr, &fdmt->states[s], false);
	fdmt->t_init.stop();

#ifdef DUMP_STATE
	// dump init state to disk
	char buf[128];
	array4d_copy_to_host(&fdmt->states[s]);
	sprintf(buf, "state_s%d.dat", 0);
	array4d_dump(&fdmt->states[s], buf);
	sprintf(buf, "initstate_e%d.dat", fdmt->execute_count);
	array4d_dump(&fdmt->states[s], buf);
#endif

	// actually execute the iterations on the GPU
	fdmt_execute_iterations(fdmt, nbeams);
	fdmt->t_update_ostate.start();
	fdmt_update_ostate(fdmt, ibeam, nbeams);
	fdmt->t_update_ostate.stop();

#ifdef DUMP_STATE
	array4d_t* currstate = &fdmt->states[fdmt->curr_state_idx];
	sprintf(buf, "finalstate_e%d.dat", fdmt->execute_count);
	array4d_copy_to_host(currstate);
	array4d_dump(currstate, buf);
	sprintf(buf, "ostate_e%d.dat", fdmt->execute_count);
	array4d_copy_to_host(&fdmt->ostate);
	array4d_dump(&fdmt->ostate, buf);
#endif
	return 0;
}

int fdmt_execute(fdmt_t* fdmt, fdmt_dtype* indata, fdmt_dtype* outdata)
{


	for(int ibeam = 0; ibeam < fdmt->nbeams; ibeam += fdmt->nbeams_alloc) {
		// do everything starting at ibeam and nbeams number
		// there must be a more elegant way to do this
		int nbeams_remaining = fdmt->nbeams - ibeam;
		assert(nbeams_remaining > 0);
		int nbeams = min(fdmt->nbeams_alloc, nbeams_remaining);
		assert(nbeams > 0);
		assert(nbeams <= fdmt->nbeams_alloc);
		fdmt_execute_batch(fdmt, indata, outdata, ibeam, nbeams);
	}

	if (fdmt->dump_data) {
		array4d_t outarray;
		outarray.d = outdata;
		outarray.nw = 1;
		outarray.nx = fdmt->nbeams;
		outarray.ny = fdmt->max_dt;
		outarray.nz = fdmt->nt;
		fdmt->t_copy_back.start();
		fdmt_copy_valid_ostate2(fdmt, &outarray);
		fdmt->t_copy_back.stop();
	}

	fdmt->execute_count += 1;
	//gpuErrchk(cudaDeviceSynchronize());

	return 0;
}

__global__ void fdmt_set_weights_kernel(const __restrict__ fdmt_dtype* ostate, fdmt_dtype* weights, int max_dt, int stride)
{
	for (int idx = blockIdx.x * blockDim.x + threadIdx.x;
			idx < max_dt;
			idx += blockDim.x * gridDim.x)
	{
		fdmt_dtype nhits = ostate[stride * idx];
	    weights[idx] = rsqrtf(nhits);
		//weights[idx] = 1.;
	    //weights[idx] = nhits;
	}
}


// Run ones through the FDMT so we can work out the weighting function
// Use the supplied inarra as working memory
int fdmt_calculate_weights(fdmt_t* fdmt)
{
	int nruns = fdmt->max_dt/fdmt->nt + 1; // add 1 for extraq giggles
	array4d_t inarr_dummy;
	inarr_dummy.nw = 1;
	inarr_dummy.nx = 1;
	inarr_dummy.ny = 1;
	inarr_dummy.nz = 1;

	for (int ii = 0; ii < nruns; ++ii) {
		int s = 0;
		fdmt_initialise_gpu(fdmt, &inarr_dummy, &fdmt->states[s], true); // Initialise state with ones
		fdmt_execute_iterations(fdmt, 1); //  actually execute the iterations on the GPU
		fdmt_update_ostate(fdmt, 0, 1); // update the output state
	}

	// set weights array
	int blocksize = 256;
	int nblocks = (fdmt->max_dt + blocksize - 1) / fdmt->max_dt;
	fdmt_set_weights_kernel<<<nblocks, blocksize>>>(fdmt->ostate.d_device, fdmt->weights.d_device, fdmt->max_dt, fdmt->max_dt + fdmt->nt);
	//gpuErrchk(cudaDeviceSynchronize());
	if (fdmt->dump_data) {
		array4d_copy_to_host(&fdmt->weights);
		array4d_dump(&fdmt->weights, "fdmtweights.dat");
	}

	// clear state
	for (int s = 0; s < 2; ++s) {
		array4d_zero(&fdmt->states[s]);
	}

	// clear ostate
	array4d_zero(&fdmt->ostate);
	//gpuErrchk(cudaDeviceSynchronize());

}

void fdmt_print_timing(fdmt_t* fdmt)
{
	cout << "FDMT Timings:" << endl;
	cout << "Initialisation: " << endl <<  fdmt->t_init << endl;
	cout << "Copy in:" << endl << fdmt->t_copy_in << endl;
	cout << "Execute: " << endl << fdmt->t_iterations << endl;
	cout << "Update Ostate: " << endl << fdmt->t_update_ostate << endl;
	cout << "Copy back: " << endl << fdmt->t_copy_back << endl;

}
