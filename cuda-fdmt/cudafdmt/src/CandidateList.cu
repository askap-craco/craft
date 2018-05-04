/*
 * CandidateList.cpp
 *
 *  Created on: 10 Feb 2017
 *      Author: ban115
 */

#include "CandidateList.h"
#include "cuda_utils.h"
#include <stdlib.h>

__global__ void printit(unsigned int *t) {
	printf("It was %d\n", *t);
}

CandidateList::CandidateList(unsigned int max_cand) {
	// Works using unified memory
	gpuErrchk(cudaMallocHost(&m_ncand, sizeof(unsigned int)));
	gpuErrchk(cudaMallocHost(&m_max_cand, sizeof(unsigned int)));
	gpuErrchk(cudaMallocHost(&m_candidates, max_cand * sizeof(candidate_t)));
	*m_max_cand = max_cand;
	clear();
}

//CandidateList::CandidateList(const& CandidateList other) :
//		m_max_cand(other.m_max_cand),
//		m_ncand(other.m_ncand),
//		m_candidates(other.m_candidates)
//{
//
//}


CandidateList::~CandidateList() {
	gpuErrchk(cudaFreeHost(m_candidates));
	gpuErrchk(cudaFreeHost(m_ncand));
	gpuErrchk(cudaFreeHost(m_max_cand));
}

void __host__ CandidateList::clear() {
	// gpuErrchk(cudaMemset(m_ncand, 0, sizeof(int)));
	*m_ncand = 0;
}

unsigned int __host__ CandidateList::ncand() {
	return *m_ncand;
}

__host__ unsigned int CandidateList::copy_to_sink(CandidateSink& sink, size_t sampno) {
	unsigned int ncand = CandidateList::ncand();
	if (ncand >= *m_max_cand-1) {
		//printf("FREDDA WARNING: Overflowed candidate buffer in block starting at sampno %d\n", sampno);
	}
	for (unsigned int i = 0; i < ncand; ++i) {
		candidate_t* c = &m_candidates[i];
		sink.add_candidate(c->ibeam, c->idt, sampno+ c->t, c->ibc, c->sn);
	}
	sink.flush();

	return ncand;
}

__device__ unsigned int CandidateList::add_candidate(int ibeam, int idt, int t, int ibc,
		float sn) {
	candidate_t* c = m_candidates + *m_ncand;
	c->ibeam = ibeam;
	c->idt = idt;
	c->t = t;
	c->ibc = ibc;
	c->sn = sn;
	return atomicInc(m_ncand, *m_max_cand);
}
