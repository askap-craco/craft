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

const size_t MAX_HDR_SIZE = 4096;

/**
 * Reads a sigproc file.
 * Sigproc order is TIF order -
 * T = time, I = IF, F= channel. F fastest
 *
 * @see http://sigproc.sourceforge.net/sigproc.pdf
 */
class SigprocFile {
public:
	SigprocFile(const char* filename);
	virtual ~SigprocFile();
	const char* header_find(const char* hname) const;
	int header_int(const char* hname) const;
	double header_double(const char* hname) const;
	size_t read_samples_uint8(size_t nt, uint8_t* output);
	double last_sample_elapsed_seconds();
	double last_sample_mjd();
	size_t seek_sample(size_t t);
	float dm_of_idt(int idt);

	double m_fch1;
	double m_foff;
	double m_tstart;
	double m_tsamp;
	int m_nifs;
	int m_nchans;
	int m_nbits;

	size_t m_samples_read;

	FILE* m_file;
	char* m_filename;
	char m_hdr[MAX_HDR_SIZE];
	size_t m_hdr_nbytes;


};

#endif /* SIGPROCFILE_H_ */
