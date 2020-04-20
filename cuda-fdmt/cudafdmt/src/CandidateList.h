/*
 * CandidateList.h
 *
 *  Created on: 10 Feb 2017
 *      Author: ban115
 */

#ifndef CANDIDATELIST_H_
#define CANDIDATELIST_H_
#include <cuda.h>

#include "CandidateSink.h"

typedef struct _candidate
{
	 int ibeam;
	 int idt;
	 int t;
	 int ibc;
	float sn;
} candidate_t;

class CandidateList {
public:
	CandidateList(unsigned int max_cand);
	//CandidateList(&CandidateList other);

	virtual ~CandidateList();

	unsigned int* m_max_cand;
	unsigned int* m_ncand;
	candidate_t* m_candidates;
	__host__ void clear();
	__host__ unsigned int ncand();
	__host__ unsigned int copy_to_host();
	__host__ unsigned int copy_to_sink(CandidateSink& sink, size_t sampno);
	__device__ unsigned int add_candidate( int ibeam,  int idt,  int t,  int ibc, float sn);
};

#endif /* CANDIDATELIST_H_ */
