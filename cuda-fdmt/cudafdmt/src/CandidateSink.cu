/*
 * CandidateSink.cpp
 *
 *  Created on: 28 Oct 2016
 *      Author: ban115
 */

#include "CandidateSink.h"
#include "SigprocFileSet.h"
#include <stdio.h>
#include <stdlib.h>

CandidateSink::CandidateSink(SigprocFileSet* srcfile, const char* filename) {
	m_srcfile = srcfile;
	if (filename == NULL) {
		char fname[4096];
		sprintf(fname, "%s.cand", srcfile->name());
		m_candfile = fopen(fname, "w");
	} else {
		m_candfile = fopen(filename, "w");
	}
	if (!m_candfile) {
		perror("Could not open candidate file");
		exit(EXIT_FAILURE);
	}
	fprintf(m_candfile, "# S/N, sampno, time from file start, boxcar, dm, beamno\n");
}

CandidateSink::~CandidateSink() {
	fflush(m_candfile);
	fclose(m_candfile);
}

void CandidateSink::add_candidate(int ibeam, int idt, int t, int ibc, float sn)
{
	// Typical line: 11.1591 7505890 9499.64 12      178     736.883 6       7505720 7507042
	// S/N sample_number dunno width dunno dm dunno dunno dunno
	// dt = 4.15ms * DM * (nu1**-2 - nu2**-2)
	// DM =

	float  dm = m_srcfile->dm_of_idt(idt);
	size_t sampno = t + m_srcfile->curr_sample();
	double time_from_file = sampno*m_srcfile->tsamp();
	fprintf(m_candfile, "%f %lu %f %d %0.3f %d\n", sn, sampno, time_from_file,
			ibc, dm, ibeam);
}
