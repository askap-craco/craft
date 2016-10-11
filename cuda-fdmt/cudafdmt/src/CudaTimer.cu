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
}

CudaTimer::~CudaTimer() {
	gpuErrchk(cudaEventDestroy(m_start));
	gpuErrchk(cudaEventDestroy(m_stop));
}

void CudaTimer::start() {
	gpuErrchk(cudaEventRecord(m_start, m_stream));
}

void CudaTimer::stop() {
	gpuErrchk(cudaEventRecord(m_stop, m_stream));
	sync_stop();

}

void CudaTimer::sync_start() {
	gpuErrchk(cudaEventSynchronize(m_start));
}

void CudaTimer::sync_stop() {
	gpuErrchk(cudaEventSynchronize(m_stop));
}

float CudaTimer::get_elapsed_time() {
	float ms;
	gpuErrchk(cudaEventElapsedTime(&ms, m_start, m_stop));
	return ms;
}
