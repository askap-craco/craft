/*
 * FileGroup.cpp
 *
 *  Created on: 1 Nov 2016
 *      Author: ban115
 */

#include <assert.h>
#include "SigprocFileSet.h"
#include "SigprocFile.h"


SigprocFileSet::SigprocFileSet(int argc, char* filenames[]) {
	m_nbeams = 0;
	for(int i = 0; i < argc; i++) {
		SigprocFile* curr_file = new SigprocFile(filenames[i]);
		m_files.push_back(curr_file);
		if (i == 0) {
			first_file = curr_file;
		} else {
			assert(first_file->nchans() == curr_file->nchans());
			assert(first_file->nbits() == curr_file->nbits());
			assert(first_file->tstart() == curr_file->tstart());
			assert(first_file->tsamp() == curr_file->tsamp());
			assert(first_file->fch1() == curr_file->fch1());
			assert(first_file->foff() == curr_file->foff());
		}
		assert(curr_file->nifs() == 1); // Otherwise the ordering will be funny
		m_nbeams += curr_file->nbeams();
	}

}

SigprocFileSet::~SigprocFileSet() {
	for(int i = 0; i < m_files.size(); i++) {
		delete m_files.at(i);
	}
}

size_t SigprocFileSet::read_samples_uint8(size_t nt, uint8_t* output)
{
	// Returns BFT ordering
	read_timer.start();
	int beamno = 0;
	size_t nread;
	for(int i = 0; i < m_files.size(); i++) {
		SigprocFile* fin = m_files.at(i);
		nread = fin->read_samples_uint8(nt, &output[beamno*nchans()*nt]);
		if (nread != nt) {
			break;
		}
		beamno += fin->nbeams();
	}
	read_timer.stop();
	return nread;
}
