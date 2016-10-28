/*
 * Boxcar functions
 */

#include "fdmt_utils.h"
#include "array.h"
#include "CandidateSink.h"

const int NBOX = 32; // Needs to be the warp size actually

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

int boxcar_do_cpu(const array4d_t* indata, array4d_t* outdata)
{
	// Inshape: [nbeams, 1, ndt, nt]
	// outshape: [nbeams, ndt, nt, nbox=32]

	int nbeams = indata->nw;
	assert(indata->nx == 1);
	int ndt = indata->ny;
	int nt = indata->nz;
	outdata->nw = nbeams;
	outdata->nx = ndt;
	outdata->ny = nt;
	outdata->nz = NBOX;
	fdmt_dtype* inp = indata->d;
	fdmt_dtype* outp = outdata->d;

	fdmt_dtype state[NBOX];
	fdmt_dtype history[NBOX];

	bzero(state, sizeof(fdmt_dtype)*NBOX);
	bzero(history, sizeof(fdmt_dtype)*NBOX);

	for(int b = 0; b < nbeams; ++b) {
		for(int idt = 0; idt < ndt; ++idt) {
			for(int t = 0; t < nt; ++t) {
				int inidx = array4d_idx(indata, b, 0, idt, t);
				fdmt_dtype v = inp[inidx];
				for (int ibc = 0; ibc < NBOX; ++ibc) {
					int history_index = (t - ibc - 1) % NBOX;
					state[ibc] += v - history[history_index];
					int outidx = array4d_idx(outdata, b, idt, t, ibc);
					outp[outidx] = state[ibc];
				}
				history[(t - 1) % NBOX] = v;
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

int boxcar_threshonly(const array4d_t* indata, fdmt_dtype thresh,
		CandidateSink& sink) {
	int nbeams = indata->nw;
	assert(indata->nx == 1);
	int ndt = indata->ny;
	int nt = indata->nz;
	int ncand = 0;

	for(int b = 0; b < nbeams; ++b) {
		for(int idt = 0; idt < ndt; ++idt) {
			for(int t = 0; t < nt; ++t) {
				int inidx = array4d_idx(indata, b, 0, idt, t);
				fdmt_dtype v = indata->d[inidx];
				if (v > thresh) {
					sink.add_candidate(b, idt, t, 0, v);
					ncand += 1;
				}
			}
		}
	}

	return ncand;
}
