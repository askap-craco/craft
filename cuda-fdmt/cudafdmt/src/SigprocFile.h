/*
 * SigprocFile.h
 *
 *  Created on: 27 Oct 2016
 *      Author: ban115
 */

#ifndef SIGPROCFILE_H_
#define SIGPROCFILE_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "DataSource.h"

const size_t MAX_HDR_SIZE = 4096;

/**
 * Reads a sigproc file.
 * Sigproc order is TIF order -
 * T = time, I = IF, F= channel. F fastest
 *
 * @see http://sigproc.sourceforge.net/sigproc.pdf
 */
class SigprocFile : public DataSource {
public:
	SigprocFile(const char* filename);
	virtual ~SigprocFile();
	const char* header_find(const char* hname) const;
	int header_int(const char* hname) const;
	double header_double(const char* hname) const;
	size_t read_samples_uint8(size_t nt, uint8_t* output);
	size_t read_samples(void** output);
	double last_sample_elapsed_seconds();
	double last_sample_mjd();

	int nifs() {
		return m_nifs;
	}
	int nbits() {
		return m_nbits;
	}
	int nbeams() {
		return nifs();
	}
	int nchans() {
		return m_nchans;
	}
	int npols() {
		return 1;
	}
	int nants() {
		return 1;
	}
	size_t samples_read() {
		return m_samples_read;
	}
	size_t curr_sample() {
		return m_curr_sample;
	}
	double fch1() {
		return m_fch1;
	}
	double foff() {
		return m_foff;
	}
	double tstart() {
		return m_tstart;
	}
	double tsamp() {
		return m_tsamp;
	}
	DataOrder data_order() {
		return DataOrder::TFBP;
	}

	char* name() {
		return m_filename;
	}

	void advise_block(off_t nt);

	size_t seek_sample(size_t t);


private:
	double m_fch1;
	double m_foff;
	double m_tstart;
	double m_tsamp;
	int m_nifs;
	int m_nchans;
	int m_nbits;

	size_t m_samples_read;
	size_t m_curr_sample;

	FILE* m_file;
	int m_fd;
	char* m_filename;
	char m_hdr[MAX_HDR_SIZE];
	size_t m_hdr_nbytes;

};

#endif /* SIGPROCFILE_H_ */
