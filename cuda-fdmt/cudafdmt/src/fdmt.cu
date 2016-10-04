
#include "cpu_kernels.h"
#include "gpu_kernels.h"
#include "cuda_fdmt.h"
#include "fdmt.h"


float dm_delay(const float f1, const float f2) {
  return 4.14e9*(isquaref(f1) - isquaref(f2));
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
  //  printf("rff %f %f %f %f %f\n", f1_start, f1_end, f2_start, f2_end, rf);

  return rf;
}

__host__ __device__ int calc_delta_t(const fdmt_t* fdmt, float f_start, float f_end)
{
  float rf = cff(f_start, f_end, fdmt->fmin, fdmt->fmax);
  float delta_tf = ((float)fdmt->max_dt-1.0) * rf;
  int delta_t = (int)ceilf(delta_tf);

  //  printf("delta t: rf %f delta_tf %f delta_t %d\n", rf, delta_tf, delta_t);
  return delta_t;
}

int fdmt_create(fdmt_t* fdmt, float fmin, float fmax, int nf, int max_dt, int nbeams)
{
  fdmt->max_dt = max_dt;
  fdmt->fmin = fmin;
  fdmt->fmax = fmax;
  fdmt->nf = nf;
  fdmt->df = (fdmt->fmax - fdmt->fmin)/((float) fdmt->nf);
  fdmt->order = (int)ceil(log(fdmt->nf)/log(2.0));
  fdmt->nbeams = nbeams;
  assert(nf > 0);
  assert(max_dt > 0);
  assert(1<<fdmt->order >= fdmt->nf);
  assert(nbeams >= 1);

  // TODO: CHeck it's important that fmin < fmax??
  assert(fmin > 0);
  assert(fmax > 0);

  //deltaT = int(np.ceil((maxDT-1) *(1./f_min**2 - 1./(f_min + deltaF)**2) / (1./f_min**2 - 1./f_max**2)))

  //fdmt->delta_t = (int)(ceilf((fdmt->maxDT-1) *(isquaref(fdmt->f_min) - isquaref(fdmt->f_min + fdmt->delta_f)) / (isquaref(f_min) - isquaref(f_max))));

  // Delta_t here is the number of time samples the maximum DM trajectory traverses
  // In the lowest channel. It is equivalent to the number of Dm trials you need to do
  // In the lowest channel to get out to the highest DM we asked for.
  fdmt->delta_t = calc_delta_t(fdmt, fdmt->fmin, fdmt->fmin + fdmt->df);
  fdmt->delta_t += 1; // Siglhtly different definition to origiinal

  // Allocate states as ping-pong buffer
  fdmt->state_size = fdmt->nbeams * fdmt->nf*fdmt->delta_t * fdmt->max_dt;
  fdmt->state_nbytes = fdmt->state_size * sizeof(fdmt_dtype);
  for (int s = 0; s < 2; s++) {
    fdmt->states[s].nw = fdmt->nbeams;
    fdmt->states[s].nx = fdmt->max_dt;
    fdmt->states[s].ny = fdmt->delta_t;
    fdmt->states[s].nz = fdmt->nf;
    array4d_malloc(&fdmt->states[s]);
  }

  return 0;
}


int fdmt_initialise(const fdmt_t* fdmt, array3d_t* indata, array4d_t* state)
{

  // indata is 3D array: (nbeams, nf, nt)
  // State is a 4D array: (nbeams, nf, deltat, nt) ( for the moment)

  assert(indata->nx == fdmt->nbeams);
  assert(indata->ny == fdmt->nf);
  assert(indata->nz == fdmt->max_dt);

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
		for (int t = 0; t < fdmt->max_dt; t++) {
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
		  int imidx = array3d_idx(indata, beam, c, indata->ny - 1);

		  // The state for dt=d = the state for dt=(d-1) + the time-reversed input sample
		  // for each time
		  // (TODO: Not including a missing overlap here)
		  for (int j = idt; j < fdmt->max_dt; j++) {
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

int fdmt_execute(fdmt_t* fdmt, fdmt_dtype* indata, fdmt_dtype* outdata)
{
  array3d_t inarr = {.nx = fdmt->nbeams, .ny = fdmt->nf, .nz = fdmt->max_dt};
  inarr.d = indata;

  // Make the final outstate - this saves a memcpy on the final iteration
  array4d_t outstate;
  outstate.nw = fdmt->states[0].nw;
  outstate.nx = fdmt->states[0].nx;
  outstate.ny = fdmt->states[0].ny;
  outstate.nz = fdmt->states[0].nz;
  outstate.d = outdata;

  // Start that puppy up
  int s = 0;
  fdmt_initialise(fdmt, &inarr, &fdmt->states[s]);
  array4d_copy_to_device(&fdmt->states[s]);

#ifdef DUMP_STATE
  char buf[128];
  sprintf(buf, "state_s%d.dat", 0);
  array4d_dump(&fdmt->states[s], buf);
#endif

  for (int iter = 1; iter < fdmt->order+1; iter++) {
    //printf("Iteration %d\n", iter);
    array4d_t* currstate = &fdmt->states[s];
    array4d_t* newstate;

    // If it's the last iteration, cheekily substitute the output pointer
    // to save a memcopy
    // if (iter == fdmt->order) {
    s = (s + 1) % 2;
    if (iter == fdmt->order) {
      newstate = &outstate;
      newstate->d_device = fdmt->states[s].d_device;
      printf("Setting outstate\n");
    } else {
      newstate = &fdmt->states[s];
    }
    fdmt_iteration(fdmt, iter, currstate, newstate);
    cuda_fdmt_iteration<<<fdmt->nbeams, fdmt->max_dt >>>(*fdmt, iter, *currstate, *newstate);
#ifdef DUMP_STATE
    array4d_copy_to_host(newstate);
    sprintf(buf, "state_s%d.dat", iter);
    array4d_dump(newstate, buf);
#endif



    //printf("Finisehd iteration %d\n", iter);

  }

  //printf("Returing form execute\n");

  return 0;
}
