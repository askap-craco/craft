/*
 * DadaSet.h
 *
 *  Created on: 25 May 2018
 *      Author: ban115
 */

#ifndef FILDIRSET_H_
#define FILDIRSET_H_
#include "DataSource.h"
#include "SigprocFileSet.h"
#include <assert.h>
#include <vector>
#include <glob.h>



class FilDirSet : public DataSource {
public:
	FilDirSet(int nt, int argc, char* dirnames[]);
	virtual ~FilDirSet();

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
		return m_latest_tstart;
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


private:
	std::vector<SigprocFileSet*> m_sources;
	SigprocFileSet* first_source;
	int m_nt;
	int m_current_ant;
	glob_t m_filglob;
	double m_latest_tstart;

};


#endif /* FILDIRSET_H_ */
