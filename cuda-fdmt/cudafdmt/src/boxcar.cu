/*
 * Boxcar functions
 */

#include "fdmt_utils.h"
#include "array.h"
#include "CandidateSink.h"
#include "boxcar.h"

int mod(int a, int b)
{
    int r = a % b;
    return r < 0 ? r + b : r;
}

__global__ void boxcar_do_kernel(const __restrict__ fdmt_dtype* indata,
		fdmt_dtype* __restrict__ outdata,
		int nt)
{
	__shared__ fdmt_dtype history[NBOX];
	int ibeam = blockIdx.x;
	int nbeams = gridDim.x;

	int idt = blockIdx.y;
	int max_dt = gridDim.y;

	int off = max_dt*(idt + ibeam*nbeams);
	const fdmt_dtype* iptr = indata + off;
	fdmt_dtype* optr = outdata + off;

	int ibc = threadIdx.x;
	int tidx = threadIdx.x;

	// initialise history
	// TODO: Load history from previous run. This will be overwritten with the state
	//history[ibc] = iptr[ibc];
	history[ibc] = 0;

	// Initialise state from history. This is basically a 'sum scan' in reverse
	// order. i.e. history[n] = sum_{i=n+1}^{NBOX}{history[i]}. Ideally you'd do a
	// Work efficient parallel scan (see http://http.developer.nvidia.com/GPUGems3/gpugems3_ch39.html)
	// But, given we're only doing 32 sums, it's probably overkill for what we want. Basically we'll
	// Sum in place in the history, and then set the thread states once we're done
	if (ibc == 0) {
		for(int t = NBOX-2; t >= 0; --t) {
			history[t] += history[t+1];
		}
	}

	__syncthreads();
	// setup the state for *this* thread (which sits in a register)
	fdmt_dtype state = history[ibc];

	// Need to load the history into shared memory again
	history[ibc] = iptr[ibc];

	for(int t = 0; t < nt; ++t) {
		// Should be a LDU instruction - global load across all threads
		fdmt_dtype v = iptr[t];
		int history_index = (t - ibc - 1) % NBOX;
		// the access to the history should have no bank conflicts, as each thread access a different bank
		state += v - history[history_index];

		// write input back to history
		if (ibc == 0) {
			history[history_index] = v;
		}

		// write state into output
		optr[ibc] = state;

		// increment output pointer
		optr += NBOX;

		__syncthreads();
	}

	// TODO: write history so we can do previous run

}

int boxcar_do_cpu(const array4d_t* indata, array4d_t* outdata, array4d_t* boxcar_history)
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

	//boxcar_history.nw = 1;
	//boxcar_history.nx = nbeams;
	//boxcar_history.ny = nd;
	//boxcar_history.nz = NBOX;

	fdmt_dtype* inp = indata->d;
	fdmt_dtype* outp = outdata->d;

	for(int b = 0; b < nbeams; ++b) {
		for(int idt = 0; idt < ndt; ++idt) {
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

			for(int t = 0; t < nt; ++t) {
				int inidx = array4d_idx(indata, b, 0, idt, t);
				fdmt_dtype vin = inp[inidx];
				for (int ibc = 0; ibc < NBOX; ++ibc) {
					int history_index = mod((-t + ibc),  NBOX);
					int outidx = array4d_idx(outdata, b, idt, t, ibc);
					state[ibc] += vin - history[history_index];
					outp[outidx] = state[ibc]/(sqrtf((float) (ibc + 1)));
				}
				int ohistidx = mod(-t-1, NBOX);
				history[ohistidx] = vin;
			}
		}
	}
}
int boxcar_do(array4d_t* indata, array4d_t* outdata)
{
	// Inshape: [nbeams, 1, ndt, nt]
	// outshape: [nbeams, ndt, nt, nbox=32]
	return 0;

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
