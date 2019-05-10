/*
 * CudaTimer.cpp
 *
 *  Created on: 4 Oct 2016
 *      Author: ban115
 */

#include "CudaTimer.h"
#include "cuda_utils.h"

CudaTimer::CudaTimer(cudaStream_t stream) {
	gpuErrchk(cudaEventCreate(&m_start));
	gpuErrchk(cudaEventCreate(&m_stop));

	m_stream = stream;
	m_total_time = 0;
	m_ncalls = 0;
}

CudaTimer::~CudaTimer() {
	gpuErrchk(cudaEventDestroy(m_start));
	gpuErrchk(cudaEventDestroy(m_stop));
}

void CudaTimer::start() {
	gpuErrchk(cudaEventRecord(m_start, m_stream));
	cputimer.start();
}

void CudaTimer::stop() {
	gpuErrchk(cudaEventRecord(m_stop, m_stream));
	sync_stop();
	m_ncalls += 1;
	float t = CudaTimer::get_elapsed_time();
	m_total_time += t;
	cputimer.stop();
}

void CudaTimer::sync_start() {
	gpuErrchk(cudaEventSynchronize(m_start));
}

void CudaTimer::sync_stop() {
	gpuErrchk(cudaEventSynchronize(m_stop));
}

float CudaTimer::get_elapsed_time() {
	// Returns elapsed time in milliseconds.
	float ms;
	if (m_ncalls == 0) { // if it hasn't been called, cudaEventElapsedTime fails with invalidResourceHandle
		ms = 0;
	} else {
		gpuErrchk(cudaEventElapsedTime(&ms, m_start, m_stop));
	}
	m_last = ms;
	return ms;
}

float CudaTimer::get_average_time() {
	// average elapsed time per clal in ms
	return m_total_time/(float)m_ncalls;
}
