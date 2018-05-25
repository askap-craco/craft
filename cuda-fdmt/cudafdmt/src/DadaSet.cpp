/*
 * DadaSet.cpp
 *
 *  Created on: 25 May 2018
 *      Author: ban115
 */

#include "DadaSet.h"


DadaSet::DadaSet(int nt, int argc, char* keynames[]) {

	for (int i = 0; i < argc; ++i) {
		DadaSource* curr_source = new DadaSource(nt, keynames[i], true);
		m_sources.push_back(curr_source);
		if (i == 0) {
			first_source = m_sources[0];
		} else {
			assert(first_source->nchans() == curr_source->nchans());
			assert(first_source->nbits() == curr_source->nbits());
			//assert(first_source->tstart() == curr_source->tstart());
			assert(first_source->tsamp() == curr_source->tsamp());
			assert(first_source->fch1() == curr_source->fch1());
			assert(first_source->foff() == curr_source->foff());
			assert(first_source->data_order() == curr_source->data_order());
		}
	}
	m_nt = nt;
	m_current_ant = 0;
	assert(m_nt > 0);
	sync();
}

DadaSet::~DadaSet() {
	for(int i = 0; i < m_sources.size(); ++i) {
		delete m_sources[i];
	}
	m_sources.clear();
}


void DadaSet::sync() {
	// Keep drawing samples until all the sources are synchronised - tries to detect if it will never happen
	double first_source_mjd = first_source->current_mjd();
	int64_t latest_sampnum = 0;

	// Loop through inputs - check they're misaligned by an integral block count, and find which one is furthest into the future
	for (int i = 0; i < m_sources.size(); ++i) {
		DadaSource* curr_source = m_sources.at(i);
		double mjd = curr_source->current_mjd();
		int64_t sampnum = curr_source->current_sample_relative_to(first_source_mjd);
		if (sampnum % m_nt != 0) {
			printf("Dada sources not synchronised. d1 %o mjdstart %f d2=%o mjdstart=%f sampnum=%ld \n",
					first_source->dada_key(), first_source_mjd, curr_source->dada_key(), mjd,
					sampnum);
			exit(EXIT_FAILURE);
		}
		if (sampnum > latest_sampnum) {
			latest_sampnum = sampnum;
		}
	}

	// max_sampdiff is the sample number we want to get to for all the input sources

	// Loop through inputs again, and read blocks until aligned with latest_sampnum
	void* dummy_output;
	for (int i = 0; i < m_sources.size(); ++i) {
		DadaSource* curr_source = m_sources.at(i);
		while(curr_source->current_sample_relative_to(first_source_mjd) < latest_sampnum) {
			curr_source->read_samples(&dummy_output);
		}
		assert(curr_source->current_sample_relative_to(first_source_mjd) == latest_sampnum);
	}

	// all sources should be synchronised at this point.
}

size_t DadaSet::read_samples(void** output)
{
	DadaSource* curr_source = m_sources.at(m_current_ant);
	m_current_ant = (m_current_ant + 1) % nants();
	return curr_source->read_samples(output);
}


