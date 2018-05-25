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
#include <assert.h>
#include <vector>



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

	size_t seek_sample(size_t t) {
		printf("Ignoring seek sample %d\n", t);
		size_t boff = 0;
//		for(int i = 0; i < m_sources.size(); i++) {
//			boff = m_sources.at(i)->seek_sample(t);
//		}
		return boff;
	}
	size_t read_samples(void** output);


private:
	std::vector<DadaSource*> m_sources;
	DadaSource* first_source;
	void sync();
	int m_nt;
	int m_current_ant;
};


#endif /* DADASET_H_ */
