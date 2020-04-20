/*
 * CandidateSink.h
 *
 *  Created on: 28 Oct 2016
 *      Author: ban115
 */

#ifndef CANDIDATESINK_H_
#define CANDIDATESINK_H_

#include "SigprocFileSet.h"
#include "SigprocFile.h"
#include <stdio.h>

class CandidateSink {
public:
	CandidateSink(SigprocFileSet* srcfile, const char* filename, bool negdm);
	virtual ~CandidateSink();
	void add_candidate(int beam, int idt, int t, int ibc, float sn);

private:
	SigprocFileSet* m_srcfile;
	FILE* m_candfile;
	bool m_negdm;
};

#endif /* CANDIDATESINK_H_ */
