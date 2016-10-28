/*
 * CandidateSink.h
 *
 *  Created on: 28 Oct 2016
 *      Author: ban115
 */

#ifndef CANDIDATESINK_H_
#define CANDIDATESINK_H_

#include "SigprocFile.h"
#include <stdio.h>

class CandidateSink {
public:
	CandidateSink(SigprocFile* srcfile, const char* filename);
	virtual ~CandidateSink();

	void add_candidate(int beam, int idt, int t, int ibc, float sn);

	SigprocFile* m_srcfile;
	FILE* m_candfile;
};

#endif /* CANDIDATESINK_H_ */
