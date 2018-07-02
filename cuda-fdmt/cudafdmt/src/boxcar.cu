/*
 * Boxcar functions
 */

#include "fdmt_utils.h"
#include "array.h"
#include "CandidateSink.h"
#include "CandidateList.h"
#include "boxcar.h"
#include <cuda_runtime_api.h>
#include <cub.cuh>

#define FULL_MASK 0xffffffff
#define WARP_SZ 32


// Modulous of a % b, but handles negative numbers
__host__ __device__ int mod(int a, int b)
{
	int r = a % b;
	return r < 0 ? r + b : r;
}

// Modulous of a % b, but handles negative numbers
__host__ __device__ fdmt_dtype mymax(fdmt_dtype a, fdmt_dtype b)
{
	return a > b ? a : b;
}

__inline__ __device__ int shfl_xor(int val, int mask)
{
#if (CUDART_VERSION >= 9000)
	return __shfl_xor_sync(FULL_MASK, val, mask);
#else
	return __shfl_xor(val, mask);
#endif
}


__inline__ __device__ int shfl(int val, int leader)
{
#if (CUDART_VERSION >= 9000)
	return __shfl_sync(FULL_MASK, val, leader);
#else
	return __shfl(val, leader);
#endif

}

__inline__ __device__ fdmt_dtype shfl_up(fdmt_dtype var, unsigned int delta, int width)
{
#if (CUDART_VERSION >= 9000)
	return __shfl_up_sync(FULL_MASK, var, delta, width);
#else
	return __shfl_up(var, delta, width);
#endif
}

__inline__ __device__ int ballot(int predicate)
{
#if (CUDART_VERSION >= 9000)
	return __ballot_sync(FULL_MASK, predicate);
#else
	return __ballot(predicate);
#endif
}

__inline__ __device__ int any(int predicate)
{
#if (CUDART_VERSION >= 9000)
	return __any_sync(FULL_MASK, predicate);
#else
	return __any(predicate);
#endif
}

__inline__ __device__ int all(int predicate)
{
#if (CUDART_VERSION >= 9000)
	return __all_sync(FULL_MASK, predicate);
#else
	return __all(predicate);
#endif
}




__device__ inline int lane_id(void) { return threadIdx.x % WARP_SZ; }
__device__ int warp_bcast(int v, int leader) { return __shfl(FULL_MASK, v, leader); }



// One of my favourite things ever
// From here: https://devblogs.nvidia.com/parallelforall/faster-parallel-reductions-kepler/
__inline__ __device__
int warpAllReduceSum(int val) {
	for (int mask = warpSize/2; mask > 0; mask /= 2)
	  val += shfl_xor(val, mask);
	return val;
}

__inline__ __device__
fdmt_dtype warpAllReduceMax(fdmt_dtype val) {
	for (int mask = warpSize/2; mask > 0; mask /= 2)
	  val = max(val, (fdmt_dtype)shfl_xor(val, mask));
	return val;
}

// Total hack because I"m too scared to do copy constructors:
__device__ unsigned int add_candidate(candidate_t* c, candidate_t* m_candidates, unsigned int* m_ncand, unsigned int m_max_cand) {

	// atomicInc will wrap once you go past m_max_cand.
	// This is dangerous
	// Increment the count by 1
	unsigned int old_ncand = atomicInc(m_ncand, m_max_cand);
	// set to max_ncand if it's wrapped
	if (old_ncand >= *m_ncand) {
		*m_ncand = m_max_cand; // Doesn't need to be synchronised
	}

	candidate_t* cnext = m_candidates + old_ncand;
	*cnext = *c;
//		printf("Added candidate ibeam=%d idt=%d ibc=%d t=%d sn=%f. Current ncand: %d old ncand %d max_cand %d\n",
//				 cnext->ibeam, cnext->idt,
//				cnext->ibc, cnext->t, cnext->sn, *m_ncand, old_ncand, *m_max_cand);

	return old_ncand;
}

__device__ void check(unsigned int* m_ncand) {
	__syncthreads();
	if (*m_ncand != 0) {
		printf("Argh! m_ncand is not zero! %d threadidx.x = %d threadidx.y = %d blockidx.x = %d blockidx.y = %d\n", *m_ncand, threadIdx.x, threadIdx.y, blockIdx.x, blockIdx.y);
	}
}

/* Shared memory for history */
__global__ void boxcar_do_kernel (
		const __restrict__ fdmt_dtype* indata,
		fdmt_dtype* __restrict__ outdata,
		fdmt_dtype* __restrict__ historydata,
		int nt,
		int ndt,
		fdmt_dtype threshold,
		candidate_t* m_candidates,
		unsigned int m_max_cand,
		unsigned int* m_ncand)
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
	//__syncthreads();

	// setup the state for *this* thread (which sits in a register hopefully from now on)
	fdmt_dtype state = history[ibc];

	// Re-read the  history back in from global memory
	history[ibc] = global_history[ibc];

	if (idt == 0) {
		//printf("initstate ibc=%d state=%f history[ibc]=%f\n", ibc, state, history[ibc]);
	}

	candidate_t cand;
	cand.t = -1; // t of the best detection so far. -1 for no curent detection ongoing
	cand.sn = -1; // vout for best detection so far.
	cand.ibc = ibc;
	cand.idt = idt;
	cand.ibeam = ibeam;


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
			//pri	ntf("ibc=%d hist_index = %d vin=%f state=%f vout=%f\n", ibc, history_index, vin, state, vout);
		}

		// write state into output
		if (outdata != NULL) {
			optr[ibc] = vout;
		}

		// increment output pointer
		optr += NBOX;

		// here is the fun bit. Find best detection over all times and boxcars

		// if in a detection:
		if (cand.t >= 0) {
			// Find out if the candidate ended this sample: do warp vote to find out if all boxcars are now below threshold
		  if (all(vout < threshold) || t == nt - 1) { // if all boxcars are below threshold, or we're at the end of the block
				// find maximum across all boxcars
				fdmt_dtype best_sn_for_ibc = warpAllReduceMax(cand.sn);
				// work out which ibc has the best vout - do a warp ballot of which ibc owns the best one
				int boxcar_mask = ballot(best_sn_for_ibc == cand.sn);
				int best_ibc = __ffs(boxcar_mask) - 1; // __ffs finds first set bit = lowsest ibc that had the all tiem best vout

				//				printf("End of candidate. t=%d sn=%f idt=%d bestsn %f mask=0x%x ibc=%d best_ibc=%d got best? %d m_ncand %d\n",
				//						cand.t, cand.sn, cand.idt, best_sn_for_ibc, boxcar_mask, cand.ibc, best_ibc, best_sn_for_ibc == cand.sn, *m_ncand);

				// if you're the winner, you get to write to memory. Lucky you!
				if (ibc == best_ibc) {
					add_candidate(&cand,  m_candidates,  m_ncand, m_max_cand);
				}

				// setup for next detection
				cand.sn = -1;
				cand.t = -1;
			} else { // detection on-going
				if (vout > cand.sn) { // keep best value for this boxcar
					cand.sn = vout;
					cand.t = t;
				}
			}
		} else { // not currently in a detection
			// do warp vote to see if any boxcars exceed threshold
		  if (any(vout >= threshold)) { // one of the boxcars has a detection that beats the threshold
				cand.sn = vout;
				cand.t = t;
			}
		}

	}

	// write history to global memory so we can do the next run of this kernel
	global_history[ibc] = history[ibc];


}

/* Use warp shuffle rather than shared memory */
__global__ void boxcar_do_kernel2 (
		const __restrict__ fdmt_dtype* indata,
		fdmt_dtype* __restrict__ outdata,
		fdmt_dtype* __restrict__ historydata,
		int nt,
		int ndt,
		fdmt_dtype threshold,
		candidate_t* m_candidates,
		unsigned int m_max_cand,
		unsigned int* m_ncand)
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

	// Specialize WarpScan for type fdmt_dtype
	typedef cub::WarpScan<fdmt_dtype> WarpScan;

	// Allocate WarpScan shared memory for DT_BLOCKS warps
	__shared__ typename WarpScan::TempStorage temp_storage[DT_BLOCKS];

	// Initialise state from history. This is basically an 'inclusive sum scan'
	// i.e. state[n] = sum_{i=0}^{n}{history[i]}. This next few lines does a
	// Work efficient parallel scan (see http://http.developer.nvidia.com/GPUGems3/gpugems3_ch39.html)
	// Implemented in CUB
	// Load data from global history. boxcar 0 holds the most recent sample, boxcar 1 a little older, etc, etc.
	fdmt_dtype state = global_history[ibc]; // here, state contains the history.

	// get CUB to do the inclusive sum scan using the temperator storage for this warp
	WarpScan(temp_storage[thread_dt]).InclusiveSum(state, state);

	// 'state' now contains the correct state for this boxcar

	// Read the previous sample (relevant for this boxcar width from global memory)
	fdmt_dtype vprev = global_history[ibc];

	candidate_t cand;
	cand.t = -1; // t of the best detection so far. -1 for no current detection ongoing
	cand.sn = -1; // vout for best detection so far.
	cand.ibc = ibc;
	cand.idt = idt;
	cand.ibeam = ibeam;

	const fdmt_dtype ibc_scale = sqrtf((float) (ibc + 1));
	threshold *= ibc_scale;// scale threshold, otherwise you have to do lots of processing per sample, which is wasteful.


	for(int t = 0; t < nt; ++t) {
		// Should be a LDU instruction - global load across all threads in this warp
		fdmt_dtype vin = *iptr;
		iptr++;

		// Add the curent sample and subtract the 'previous' one
		state += vin - vprev;

		// shift previous values one thread to the right. leaves vprev for ibc=0 unchanged.
		vprev = shfl_up(vprev, 1, NBOX);

		// set vprev to vin for ibc=0
		if (ibc == 0) {
			vprev = vin;
		}

		// scale output value to have constant variance per boxcar size
		//fdmt_dtype vout = state/(sqrtf((float) (ibc + 1)));
		//fdmt_dtype vout = state * ibc_scale;
		fdmt_dtype vout = state;

		// write state into output
		if (outdata != NULL) {
			optr[ibc] = vout/ibc_scale;
			// increment output pointer
			optr += NBOX;
		}

		// here is the fun bit. Find best detection over all times and boxcars

		// if in a detection:
		if (cand.t >= 0) {
			// Find out if the candidate ended this sample: do warp vote to find out if all boxcars are now below threshold
		  if (all(vout < threshold) || t == nt - 1) { // if all boxcars are below threshold, or we're at the end of the block
				// find maximum across all boxcars
				fdmt_dtype scaled_sn = cand.sn/ibc_scale; // scale by boxcar width, so all boxcars have the same variance
				fdmt_dtype best_sn_for_ibc = warpAllReduceMax(scaled_sn);
				// work out which ibc has the best vout - do a warp ballot of which ibc owns the best one
				int boxcar_mask = ballot(best_sn_for_ibc == scaled_sn);
				int best_ibc = __ffs(boxcar_mask) - 1; // __ffs finds first set bit = lowsest ibc that had the all tiem best vout

				// if you're the winner, you get to write to memory. Lucky you!
				if (ibc == best_ibc) {
					cand.sn = scaled_sn;
					// rescale sn, which is currently unrescaled
					add_candidate(&cand,  m_candidates,  m_ncand, m_max_cand);
				}

				// setup for next detection
				cand.sn = -1;
				cand.t = -1;
			} else { // detection on-going
				if (vout > cand.sn) { // keep best value for this boxcar
					cand.sn = vout;
					cand.t = t;
				}
			}
		} else { // not currently in a detection
			// do warp vote to see if any boxcars exceed threshold
		  if (any(vout >= threshold)) { // one of the boxcars has a detection that beats the threshold
				cand.sn = vout;
				cand.t = t;
			}
		}

	}

	// write vprev to global memory so we can do the next run of this kernel
	global_history[ibc] = vprev;

}

/* Use warp shuffle for history, shared memory for input */
__launch_bounds__(32, 32)
__global__ void boxcar_do_kernel3 (
		const __restrict__ fdmt_dtype* indata,
		fdmt_dtype* __restrict__ outdata,
		fdmt_dtype* __restrict__ historydata,
		fdmt_dtype* __restrict__ discarddata,
		int nt,
		int ndt,
		fdmt_dtype threshold,
		candidate_t* m_candidates,
		unsigned int m_max_cand,
		unsigned int* m_ncand,
		int maxbc)
{
	// gridDim.x = nbeams
	// gridDim.y = ndt_blocks // number of blocks of dispersion trials
	// blockDim.x = NBOX
	// blockDim.y = DT_BLOCKS = number of dispersion trials to do at once
	// TOTAL nthreads = NBOX * DT_BLOCKS
	// Assume ndt is an integer multiple of DT_BLOCKS
	// indata.shape = [nbeams, 1, ndt, ndt + nt]
	// outdata.shape = [nbeams,ndt, nt, nbox]
	// history.shape = [1, nbeams, ndt, NBOX]
	// discard data.shape [1,1, nbeams,ndt]

	int nbeams = gridDim.x;

	int ibeam = blockIdx.x; // beam index
	int grid_dt = blockIdx.y; // delta_t of grid
	int thread_dt = threadIdx.y; // delta_t from thread
	int idt = thread_dt + DT_BLOCKS*grid_dt; // total delta_t
	int ibc = threadIdx.x; // boxcar width in samples. ibc=0 is 1 sample wide.

	int in_off = array4d_idx(nbeams, 1, ndt, ndt + nt, ibeam, 0, idt, 0); // input offset
	int out_off = array4d_idx(nbeams, ndt, nt, NBOX, ibeam, idt, 0, 0); // output offset
	int hist_off = array4d_idx(1, nbeams, ndt, NBOX, 0, ibeam, idt, 0); // history offset
	int discard_off = array4d_idx(1,1,nbeams,ndt,0,0,ibeam,idt); // discard offset

	// Shared memory
	__shared__ fdmt_dtype thread_indata[DT_BLOCKS][NBOX];

	const fdmt_dtype* iptr = indata + in_off + ibc;
	fdmt_dtype* global_history = historydata + hist_off;
	fdmt_dtype* optr = outdata + out_off;

	// Specialize WarpScan for type fdmt_dtype
	typedef cub::WarpScan<fdmt_dtype> WarpScan;

	// Allocate WarpScan shared memory for DT_BLOCKS warps
	__shared__ typename WarpScan::TempStorage temp_storage[DT_BLOCKS];

	// Initialise state from history. This is basically an 'inclusive sum scan'
	// i.e. state[n] = sum_{i=0}^{n}{history[i]}. This next few lines does a
	// Work efficient parallel scan (see http://http.developer.nvidia.com/GPUGems3/gpugems3_ch39.html)
	// Implemented in CUB
	// Load data from global history. boxcar 0 holds the most recent sample, boxcar 1 a little older, etc, etc.
	fdmt_dtype state = global_history[ibc]; // here, state contains the history.

	// get CUB to do the inclusive sum scan using the temperator storage for this warp
	WarpScan(temp_storage[thread_dt]).InclusiveSum(state, state);

	// 'state' now contains the correct state for this boxcar

	// Read the previous sample (relevant for this boxcar width from global memory)
	fdmt_dtype vprev = global_history[ibc];

	candidate_t cand;
	cand.t = -1; // t of the best detection so far. -1 for no current detection ongoing
	cand.sn = -1; // vout for best detection so far.
	cand.ibc = ibc;
	cand.idt = idt;
	cand.ibeam = ibeam;

	const fdmt_dtype ibc_scale = sqrtf((float) (ibc + 1));
	threshold *= ibc_scale;// scale threshold, otherwise you have to do lots of processing per sample, which is wasteful.
	int fullt = 0;

	for(int blk = 0; blk < nt/NBOX; ++blk) {
		// Should be a coalesced
		fdmt_dtype vin = *iptr;
		iptr += NBOX;

		// store in shared memory - should be coalesced
		thread_indata[thread_dt][ibc] = vin;
		__syncthreads(); // actually should do a __syncwarp here - cuda9 required

		// loop through shared memory
		for (int t = 0 ; t < NBOX; ++t) {
			vin = thread_indata[thread_dt][t]; // 'broadcast' from shared memory
			// Add the current sample and subtract the 'previous' one
			state += vin - vprev;

			// shift previous values one thread to the right. leaves vprev for ibc=0 unchanged.
			vprev = shfl_up(vprev, 1, NBOX);

			// set vprev to vin for ibc=0
			if (ibc == 0) {
				vprev = vin;
			}

			// scale output value to have constant variance per boxcar size
			//fdmt_dtype vout = state/(sqrtf((float) (ibc + 1)));
			//fdmt_dtype vout = state * ibc_scale;
			fdmt_dtype vout = state;

			// write state into output
			if (outdata != NULL) {
				optr[ibc] = vout/ibc_scale;
				// increment output pointer
				optr += NBOX;
			}


			// here is the fun bit. Find best detection over all times and boxcars
			// if in a detection:
			if (cand.t >= 0) {

				if (vout > cand.sn) { // keep best value for this boxcar. Put here in case the final sample (next if statement) is true
					cand.sn = vout;
					cand.t = fullt;
				}
				// Find out if the candidate ended this sample: do warp vote to find out if all boxcars are now below threshold
				if (all(vout < threshold) || fullt == nt - 1) { // if all boxcars are below threshold, or we're at the end of the block
					// find maximum across all boxcars
					fdmt_dtype scaled_sn = cand.sn/ibc_scale; // scale by boxcar width, so all boxcars have the same variance
					fdmt_dtype best_sn_for_ibc = warpAllReduceMax(scaled_sn);
					// work out which ibc has the best vout - do a warp ballot of which ibc owns the best one
					int boxcar_mask = ballot(best_sn_for_ibc == scaled_sn);
					int best_ibc = __ffs(boxcar_mask) - 1; // __ffs finds first set bit = lowsest ibc that had the all tiem best vout

					// if you're the winner, you get to write to memory. Lucky you!
					if (ibc == best_ibc) {
						if (ibc <= maxbc) {
							cand.sn = scaled_sn	;
							// rescale sn, which is currently unrescaled
							add_candidate(&cand,  m_candidates,  m_ncand, m_max_cand);
						} else {
							// Global read/write by a single thread - bleah, but oh well.
							// todo: for speed keep in a register and do a warp sum at the end....
							discarddata[discard_off] += 1.0f;
						}
					}

					// setup for next detection
					cand.sn = -1;
					cand.t = -1;
				} else { // detection on-going

				}
			} else { // not currently in a detection
				// do warp vote to see if any boxcars exceed threshold
			  if (any(vout >= threshold)) { // one of the boxcars has a detection that beats the threshold
					cand.sn = vout;
					cand.t = fullt;
				}
			}

			fullt++;
		}
	}

	// write vprev to global memory so we can do the next run of this kernel
	global_history[ibc] = vprev;

}

__global__ void test_cub_prefix_sum() {

	// Specialize WarpScan for type fdmt_dtype
	typedef cub::WarpScan<fdmt_dtype> WarpScan;

	// Allocate WarpScan shared memory for one warp
	__shared__ typename WarpScan::TempStorage temp_storage;

	// Only the first warp performs a prefix sum
	if (threadIdx.x < 32)
	{
		// Obtain one input item per thread
		fdmt_dtype thread_data = (fdmt_dtype)threadIdx.x;
		// Compute warp-wide prefix sums
		printf("Before threadidx.x %d data %f\n", threadIdx.x, thread_data);
		WarpScan(temp_storage).InclusiveSum(thread_data, thread_data);
		printf("After threadidx.x %d data %f\n", threadIdx.x, thread_data);

		thread_data = (fdmt_dtype)threadIdx.x;
		// Compute warp-wide prefix sums
		WarpScan(temp_storage).ExclusiveSum(thread_data, thread_data);
		printf("After ex  threadidx.x %d data %f\n", threadIdx.x, thread_data);
	}
}


int boxcar_do_gpu(const array4d_t* indata,
		array4d_t* boxcar_data,
		array4d_t* boxcar_history,
		array4d_t* boxcar_discards,
		fdmt_dtype thresh, int max_ncand_per_block, int mindm, int maxbc,
		CandidateList* sink)
{
	// indata is the FDMT ostate: i.e. Inshape: [nbeams, 1, ndt, ndt]
	// boxcar_data shape: [nbeams, ndt, nt, nbox=32] - device and host pointers null if we don't want to save to memory
	// But we might only want to boxcar the first nt of it
	int nbeams = indata->nw;
	int nt = boxcar_data->ny;
	int ndt = indata->ny;
	assert(indata->nw  == nbeams);
	assert(indata->nx == 1);
	assert(indata->ny == ndt);
	assert(indata->nz == ndt + nt);
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
	assert(boxcar_history->d_device != NULL);
	sink->clear();
	assert(sink->ncand() == 0);

	boxcar_do_kernel3<<<grid_shape, block_shape>>>(
			indata->d_device,
			boxcar_data->d_device,
			boxcar_history->d_device,
			boxcar_discards->d_device,
			nt, ndt, thresh,
			sink->m_candidates,
			*sink->m_max_cand,
			sink->m_ncand,
			maxbc);

	gpuErrchk(cudaDeviceSynchronize());
	return 0;

}

int boxcar_do_cpu(const array4d_t* indata,
		array4d_t* outdata,
		array4d_t* boxcar_history,
		fdmt_dtype thresh, size_t sampno,
		int max_ncand_per_block, int mindm, int maxbc,
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
