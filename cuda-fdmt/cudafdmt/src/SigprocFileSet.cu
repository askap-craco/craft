/*
 * FileGroup.cpp
 *
 *  Created on: 1 Nov 2016
 *      Author: ban115
 */

#include <assert.h>
#include "SigprocFileSet.h"
#include "SigprocFile.h"


SigprocFileSet::SigprocFileSet(int nt, int argc, char* filenames[]) : m_nt(nt) {
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

	// Create read buffer
	int in_num_elements = nt*nchans()*nbeams(); // Number of elements per block read
	size_t in_num_bytes = sizeof(uint8_t)*8*in_num_elements/nbits(); //
	uint8_t* read_buf = (uint8_t*) malloc(in_num_bytes);
	assert(read_buf);

}

SigprocFileSet::~SigprocFileSet() {
	for(int i = 0; i < m_files.size(); i++) {
		delete m_files.at(i);
	}
	if (read_buf) {
		free(read_buf);
	}
}

size_t SigprocFileSet::read_samples(void** output)
{
	m_read_timer.start();
	// Returns BFT ordering
	int beamno = 0;
	size_t nread;
	size_t bytes_per_block = nchans()*m_nt*8/nbits();
	for(int i = 0; i < m_files.size(); i++) {
		SigprocFile* fin = m_files.at(i);
		nread = fin->read_samples_uint8(m_nt, &read_buf[beamno*bytes_per_block]);
		if (nread != bytes_per_block) {
			break;
		}
		beamno += fin->nbeams();
	}
	m_read_timer.stop();
	*output = read_buf;
	return nread;
}
