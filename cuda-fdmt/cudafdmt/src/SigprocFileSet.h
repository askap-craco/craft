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
#include "CpuTimer.h"


//class FileGroup : public fdmt::DataSource { // GOD I hate C++ classes
class SigprocFileSet {

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

	size_t seek_sample(int t) {
		size_t boff;
		for(int i = 0; i < m_files.size(); i++) {
			boff = m_files.at(i)->seek_sample(t);
		}
		return boff;
	}

	size_t read_samples_uint8(size_t nt, uint8_t* output);

	CpuTimer read_timer;


private:
	std::vector<SigprocFile*> m_files;
	int m_nbeams;
	SigprocFile* first_file;

};

#endif /* FILEGROUP_H_ */
