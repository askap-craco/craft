/*
 * CandidateSink.cpp
 *
 *  Created on: 28 Oct 2016
 *      Author: ban115
 */

#include "CandidateSink.h"
#include <stdio.h>
#include <stdlib.h>

CandidateSink::CandidateSink(SigprocFile* srcfile, const char* filename) {
	m_srcfile = srcfile;
	if (filename == NULL) {
		char fname[4096];
		sprintf(fname, "%s.cand", srcfile->m_filename);
		m_candfile = fopen(fname, "w");
	} else {
		m_candfile = fopen(filename, "w");
	}
	if (!m_candfile) {
		perror("Could not open candidate file");
		exit(EXIT_FAILURE);
	}
}

CandidateSink::~CandidateSink() {
	fclose(m_candfile);
}

void CandidateSink::add_candidate(int ibeam, int idt, int t, int ibc, float sn)
{
	// Typical thing: 11.1591 7505890 9499.64 12      178     736.883 6       7505720 7507042
	// S/N sample_number dunno width dunno dm dunno dunno dunno
	// dt = 4.15ms * DM * (nu1**-2 - nu2**-2)
	// DM =

	float  dm = m_srcfile->dm_of_idt(idt);
	fprintf(m_candfile, "%f %d %f %d %d %0.3f %d %d %d %d\n", sn, t + m_srcfile->m_samples_read,
				0, ibc, 0, dm, 0, 0, 0, ibeam);
}
