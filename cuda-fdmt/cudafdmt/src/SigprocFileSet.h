/*
 * FileGroup.h
 *
 *  Created on: 1 Nov 2016
 *      Author: ban115
 */

#ifndef FILEGROUP_H_
#define FILEGROUP_H_

#include <vector>
#include "DataSource.h"
#include "SigprocFile.h"


class SigprocFileSet : public DataSource {

public:
	SigprocFileSet(int argc, char* filenames[]);
	virtual ~SigprocFileSet();

	float dm_of_idt(int idt) {
		return first_file->dm_of_idt(idt);
	}
	double fch1() {
		return first_file->fch1();
	}
	double foff() {
		return first_file->foff();
	}
	char* name() {
		return first_file->name(); // TODO: Replace with siomething more meaningful
	}
	int nbeams() {
		return m_nbeams;
	}
	int nchans() {
		return first_file->nchans();
	}
	int npols() {
		return first_file->npols();
	}
	int nants() {
		return 1;
	}
	double tsamp() {
		return first_file->tsamp();
	}
	double tstart() {
		return first_file->tstart();
	}
	size_t samples_read() {
		return first_file->samples_read();
	}
	size_t curr_sample() {
		return first_file->curr_sample();
	}

	size_t seek_sample(size_t t) {
		size_t boff;
		for(int i = 0; i < m_files.size(); i++) {
			boff = m_files.at(i)->seek_sample(t);
		}
		return boff;
	}

	size_t read_samples_uint8(size_t nt, uint8_t* output);

private:
	std::vector<SigprocFile*> m_files;
	int m_nbeams;
	SigprocFile* first_file;

};

#endif /* FILEGROUP_H_ */
