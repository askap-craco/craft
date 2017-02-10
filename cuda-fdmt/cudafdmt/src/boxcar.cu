/*
 * Boxcar functions
 */

#include "fdmt_utils.h"
#include "array.h"
#include "CandidateSink.h"
#include "boxcar.h"

// Modulous of a % b, but handles negative numbers
__host__ __device__ int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}

__global__ void boxcar_do_kernel (
		const __restrict__ fdmt_dtype* indata,
		fdmt_dtype* __restrict__ outdata,
		fdmt_dtype* __restrict__ historydata,
		int nt,
		int ndt)
{
	// gridDim.x = nbeams
	// gridDim.y = ndt_blocks // number of blocks of dispersion trials
    // blockDim.x = NBOX
	// blockDim.y = DT_BLOCKS = number of dispersion trials to do at once
	// TOTAL nthreads = NBOX * DT_BLOCKS
	// Assume ndt is an integer multiple of DT_BLOCKS
	// indata.shape = [nbeams, 1, ndt, ndt]
	// outdata.shape = [nbeams,ndt, nt, nbox]
	// history.shape = [1, nbeams, ndt, NBOX]

	__shared__ fdmt_dtype thread_history[DT_BLOCKS][NBOX];
	int nbeams = gridDim.x;

	int ibeam = blockIdx.x; // beam index
	int grid_dt = blockIdx.y; // delta_t of grid
	int thread_dt = threadIdx.y; // delta_t from thread
	int idt = thread_dt + DT_BLOCKS*grid_dt; // total delta_t
	int ibc = threadIdx.x; // boxcar width in samples. ibc=0 is 1 sample wide.

	int in_off = array4d_idx(nbeams, 1, ndt, ndt, ibeam, 0, idt, 0); // input offset
	int out_off = array4d_idx(nbeams, ndt, nt, NBOX, ibeam, idt, 0, 0); // output offset
	int hist_off = array4d_idx(1, nbeams, ndt, NBOX, 0, ibeam, idt, 0); // history offset

	const fdmt_dtype* iptr = indata + in_off;
	fdmt_dtype* global_history = historydata + hist_off;
	fdmt_dtype* optr = outdata + out_off;
	fdmt_dtype* history = thread_history[thread_dt];

	// Initialise state from history. This is basically a 'sum scan' in reverse
	// order. i.e. state[n] = sum_{i=n+1}^{NBOX}{history[i]}. Ideally you'd do a
	// Work efficient parallel scan (see http://http.developer.nvidia.com/GPUGems3/gpugems3_ch39.html)
	// But, given we're only doing NBOX sums (relative to how much work we're about to do)
	// such a scan is probably overkill for what we want. Basically we'll
	// We'll borrow the history buffer for a bit to initialise the state. Or see CUB
	// history increases to the left.

	if (ibc == 0) {
		history[0] = global_history[0];

		for(int it = 1; it < NBOX; ++it) {
			history[it] = history[it - 1] + global_history[it];
			if (idt==0) {
				//printf("statecalc it=%d h[%d]=%f h[%d]=%f gh[%d]=%f off=%d\n", it, it, history[it], it-1, history[it-1], it, global_history[it], in_off);
			}
		}
	}

	// not strictly required, as all the communication is between threads of a warp, but makes me feel good.
	__syncthreads();

	// setup the state for *this* thread (which sits in a register hopefully from now on)
	fdmt_dtype state = history[ibc];

	// Re-read the  history back in from global memory
	history[ibc] = global_history[ibc];

	if (idt == 0) {
		//printf("initstate ibc=%d state=%f history[ibc]=%f\n", ibc, state, history[ibc]);
	}


	for(int t = 0; t < nt; ++t) {
		// Should be a LDU instruction - global load across all threads in this warp
		fdmt_dtype vin = iptr[t];

		// calculate which part of the history we should find our data
		int history_index = mod((-t + ibc),  NBOX);
		// the access to the history should have no bank conflicts, as each thread access a different bank
		state += vin - history[history_index];

		// write input back to history
		if (ibc == 0) {
			int ohistidx = mod(-t-1, NBOX);
			history[ohistidx] = vin;
		}

		// scale output value to have constant variance per boxcar size
		fdmt_dtype vout = state/(sqrtf((float) (ibc + 1)));

		if (idt == 0) {
			//printf("ibc=%d hist_index = %d vin=%f state=%f vout=%f\n", ibc, history_index, vin, state, vout);
		}

		// write state into output
		optr[ibc] = vout;

		// increment output pointer
		optr += NBOX;

		// you can do a __syncthreads here if you really want, but because we're wroking in a warp, and there's
		// no communication between warps (they're on different IDT indeces), we don't need to
		//__syncthreads();
	}

	// write history to global memory so we can do the next run of this kernel
	global_history[ibc] = history[ibc];


}


int boxcar_do_gpu(const array4d_t* indata,
		array4d_t* boxcar_data,
		array4d_t* boxcar_history,
		size_t sampno,
		fdmt_dtype thresh, int max_ncand_per_block, int mindm, int maxbc,
		CandidateSink& sink)
{
	// indata is the FDMT ostate: i.e. Inshape: [nbeams, 1, ndt, ndt]
	// boxcar_data shape: [nbeams, ndt, nt, nbox=32]
	// But we might only want to boxcar the first nt of it
	int nbeams = indata->nw;
	int nt = boxcar_data->ny;
	int ndt = indata->ny;
	assert(indata->nw  == nbeams);
	assert(indata->nx == 1);
	assert(indata->ny == ndt);
	assert(indata->nz == ndt);
	assert(boxcar_data->nw == nbeams);
	assert(boxcar_data->nx == ndt);
	assert(boxcar_data->ny == nt);
	assert(boxcar_data->nz == NBOX);
	assert(boxcar_history->nw == 1);
	assert(boxcar_history->nx == nbeams);
	assert(boxcar_history->ny == ndt);
	assert(boxcar_history->nz == NBOX);

	assert(mindm < ndt);

	assert(ndt % DT_BLOCKS == 0); // otherwise kernel doesn't work
	int ndt_blocks = ndt / DT_BLOCKS;
	dim3 block_shape( NBOX, DT_BLOCKS);
	dim3 grid_shape(nbeams, ndt_blocks);

	assert(indata->d_device != NULL);
	assert(boxcar_data->d_device != NULL);
	assert(boxcar_history->d_device != NULL);
	boxcar_do_kernel<<<grid_shape, block_shape>>>(
			indata->d_device,
			boxcar_data->d_device,
			boxcar_history->d_device,
			nt, ndt);
	return 0;

}

int boxcar_do_cpu(const array4d_t* indata,
		array4d_t* outdata,
		array4d_t* boxcar_history,
		size_t sampno,
		fdmt_dtype thresh, int max_ncand_per_block, int mindm, int maxbc,
		CandidateSink& sink)
{
	// Inshape: [nbeams, 1, ndt, nt]
	// outshape: [nbeams, ndt, nt, nbox=32]
	// KB checked this code on 10 Dec 2016 late at night with much pain - and it works.

	int nbeams = indata->nw;
	assert(indata->nx == 1);
	int ndt = indata->ny;
	int nt = indata->nz;
	outdata->nw = nbeams;
	outdata->nx = ndt;
	outdata->ny = nt;
	outdata->nz = NBOX;

	assert(mindm < ndt);
	//boxcar_history.nw = 1;
	//boxcar_history.nx = nbeams;
	//boxcar_history.ny = nd;
	//boxcar_history.nz = NBOX;

	fdmt_dtype* inp = indata->d;
	fdmt_dtype* outp = outdata->d;
	int ncand = 0;
	#pragma omp parallel for shared(ncand, inp, outp)
	for(int b = 0; b < nbeams; ++b) {
		int beam_ncand = 0;
		for(int idt = mindm; idt < ndt; ++idt) {
			// Break out of loop if we've exceeded ncand per block
			if (beam_ncand >= max_ncand_per_block) {
				break;
			}
			// initialise state from boxcar history
			fdmt_dtype state[NBOX];
			int histidx = array4d_idx(boxcar_history, 0, b, idt, 0);
			fdmt_dtype* history = &boxcar_history->d[histidx];

			// history increases to the left
			state[0] = history[0];
			for(int ibc = 1; ibc < NBOX; ++ibc) {
				state[ibc] = state[ibc-1] + history[ibc];
			}

			assert(state[0] == history[0]);

			// A candidate starts when any boxcar exceeds the threshold
			// and ends when all boxcars are below the threshold
			// The highest S/N, boxcar and time of the candidate are all recorded
			// and reported when the candidate finishes, or the end of the block arrives

			// Initialise candidate grouping variables
			int best_ibc = -1; // best boxcar index
			fdmt_dtype best_ibc_sn = -1; // S/N of best boxcar
			int best_ibc_t = -1; // t of best boxcar
			int cand_tstart = -1; // t when threshold first exceeded. cand_tstart >=0 signifies a candidate in progress

			for(int t = 0; t < nt; ++t) {
				int inidx = array4d_idx(indata, b, 0, idt, t);
				fdmt_dtype vin = inp[inidx];
				for (int ibc = 0; ibc < NBOX; ++ibc) {
					// Calculate boxcar value
					int history_index = mod((-t + ibc),  NBOX);
					state[ibc] += vin - history[history_index];
					// Scale boxcar value with sqrt(length) to rescale to S/N
					fdmt_dtype vout = state[ibc]/(sqrtf((float) (ibc + 1)));
					int outidx = array4d_idx(outdata, b, idt, t, ibc);
					outp[outidx] = vout;
					if (vout > thresh) { // if this  S/N is above the threshold

						//printf("boxcar exceeded threshold b/idt/t/ibc %d/%d/%d/%d vout=%f cand_tstart=%d best ibc/sn/t = %d/%f/%d\n",
								//b, idt, t, ibc, vout, cand_tstart, best_ibc, best_ibc_sn, best_ibc_t);
						// If there's no current candidate, we start one
						if (cand_tstart < 0) {
							cand_tstart = t;
						}

						// Keep the best S/N and ibc
						if (vout > best_ibc_sn) {
							best_ibc_sn = vout;
							best_ibc = ibc;
							best_ibc_t = t;
						}
					}

					// If there was a candidate,
					// AND it's ended because no boxcars exceeded the threshold this time
					// OR if it's the last sample of the block (// because I can't be bothered keeping all the best_* states)
					// THEN write out the candidate
					if ((cand_tstart >= 0) && ((best_ibc < 0) || (t == nt - 1))) {

						// Ignore this candidate if it's best boxcar is > maxbc
						if (best_ibc <= maxbc) {
#pragma omp critical
							{
								assert(best_ibc_sn > 0);
								// record candidate details (the first boxcar is 1 sample wide)
								sink.add_candidate(b, idt, best_ibc_t + sampno, best_ibc + 1, best_ibc_sn);
								++ncand;
							}
							++beam_ncand;
						}
						// reset stuff for the next candidate
						best_ibc_sn = -1;
						best_ibc = -1;
						best_ibc_t  = -1;
						cand_tstart = -1;
					}
				}


				int ohistidx = mod(-t-1, NBOX);
				history[ohistidx] = vin;
			}
		}
	}

	return ncand;
}

int boxcar_threshonly(const array4d_t* indata, size_t sampno, fdmt_dtype thresh, int max_ncand_per_block, int mindm,
		CandidateSink& sink) {
	int nbeams = indata->nw;
	assert(indata->nx == 1);
	int ndt = indata->ny;
	int nt = indata->nz;
	int ncand = 0;

#pragma omp parallel for shared(ncand)
	for(int b = 0; b < nbeams; ++b) {
		for(int idt = mindm; idt < ndt; ++idt) {
			int off = array4d_idx(indata, b, 0, idt, 0);
			for(int t = 0; t < nt; ++t) {
				int inidx = off + t;
				fdmt_dtype v = indata->d[inidx];
				if (v > thresh && ncand < max_ncand_per_block) {
#pragma omp critical
					{
						sink.add_candidate(b, idt, t+sampno, 0, v);
						ncand += 1;
					}
				}
			}
		}
		if (ncand >= max_ncand_per_block) {

		}
	}

	return ncand;
}
