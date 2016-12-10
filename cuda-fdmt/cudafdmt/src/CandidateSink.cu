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

CandidateSink::CandidateSink(SigprocFileSet* srcfile, const char* filename, bool negdm) {
	m_srcfile = srcfile;
	m_negdm = negdm;
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
	fprintf(m_candfile, "# S/N, sampno, secs from file start, boxcar, idt, dm, beamno\n");
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
	double time_from_file = t*m_srcfile->tsamp();
	if (m_negdm) {
		idt = -idt;
		dm = -dm;
	}
	fprintf(m_candfile, "%f %lu %f %d %d %0.3f %d\n", sn, t, time_from_file,
			ibc, idt, dm, ibeam);
}
