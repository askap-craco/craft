/*
 * FilDirSet.cpp
 *
 *  Created on: 25 May 2018
 *      Author: ban115
 */

#include "FilDirSet.h"
#include "InvalidSourceFormat.h"


FilDirSet::FilDirSet(int nt, int argc, char* dirnames[]) {

	for (int i = 0; i < argc; ++i) {
		char glob_pattern[1024];
		snprintf(glob_pattern, 1024, "%s/*.fil", dirnames[i]);
		glob(glob_pattern, GLOB_TILDE, NULL, &m_filglob);
		if (m_filglob.gl_pathc == 0) {
			printf("No filterbanks in %s\n", dirnames[i]);
			throw InvalidSourceFormat();
		}
		SigprocFileSet* curr_source = new SigprocFileSet(nt, m_filglob.gl_pathc, m_filglob.gl_pathv);
		m_sources.push_back(curr_source);
		if (i == 0) {
			first_source = m_sources[0];
			m_latest_tstart = first_source->tstart();
		} else {
			assert(first_source->nchans() == curr_source->nchans());
			assert(first_source->nbits() == curr_source->nbits());
			//assert(first_source->tstart() == curr_source->tstart());
			assert(first_source->tsamp() == curr_source->tsamp());
			assert(first_source->fch1() == curr_source->fch1());
			assert(first_source->foff() == curr_source->foff());
			assert(first_source->data_order() == curr_source->data_order());
			if (curr_source->tstart() > m_latest_tstart) {
				m_latest_tstart = curr_source->tstart();
			}
		}
	}

	m_nt = nt;
	m_current_ant = 0;
	assert(m_nt > 0);

	// align all sources
	seek_sample(0);

}

FilDirSet::~FilDirSet() {
	for(int i = 0; i < m_sources.size(); ++i) {
		delete m_sources[i];
	}
	m_sources.clear();
	globfree(&m_filglob);
}

size_t FilDirSet::seek_sample(size_t offset) {
	// need to offset each source to a different offset

	for(auto source : m_sources) {
		double mjddiff = m_latest_tstart - source->tstart();
		assert(mjddiff >= 0);
		int64_t sampoff = (int64_t) round(mjddiff * 86400.0 / source->tsamp());
		source->seek_sample(sampoff);
	}
	return offset;
}


size_t FilDirSet::read_samples(void** output)
{
	DataSource* curr_source = m_sources.at(m_current_ant);
	m_current_ant = (m_current_ant + 1) % nants();
	return curr_source->read_samples(output);
}

size_t FilDirSet::read_samples_ant(void** output, int iant)
{
	DataSource* curr_source = m_sources.at(iant);
	return curr_source->read_samples(output);
}


