/*
 * DadaSet.h
 *
 *  Created on: 25 May 2018
 *      Author: ban115
 */

#ifndef DADASET_H_
#define DADASET_H_
#include "DataSource.h"
#include "DadaSource.h"
#include <string>
#include <assert.h>
#include <vector>

using std::string;

class DadaSet : public DataSource {
public:
	DadaSet(int nt, int argc, char* keynames[]);
	virtual ~DadaSet();

	float dm_of_idt(int idt) {
		return first_source->dm_of_idt(idt);
	}
	double fch1() {
		return first_source->fch1();
	}
	double foff() {
		return first_source->foff();
	}
	char* name() {
		return first_source->name();
	}
	int nbits() {
		return first_source->nbits();
	}
	int nbeams() {
		return first_source->nbeams();
	}
	int nchans() {
		return first_source->nchans();
	}
	int npols() {
		return first_source->npols();
	}
	int nants() {
		return m_sources.size();
	}
	double tsamp() {
		return first_source->tsamp();
	}
	double tstart() {
		return first_source->tstart();
	}
	DataOrder data_order() {
		return first_source->data_order();
	}

	size_t current_sample() {
		return first_source->current_sample();
	}

	size_t seek_sample(size_t t);

	size_t read_samples(void** output);

	size_t read_samples_ant(void** output, int iant);
	const char* antenna_name();

	DadaSource* get_source_at(int idx) {
		return m_sources.at(idx);
	}


private:
	std::vector<DadaSource*> m_sources;
	DadaSource* first_source;
	void sync(size_t offset);
	int m_nt;
	int m_current_ant;
	std::string m_antenna_name;
};


#endif /* DADASET_H_ */
