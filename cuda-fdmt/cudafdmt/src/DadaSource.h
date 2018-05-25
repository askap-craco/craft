/*
 * DadaSource.h
 *
 *  Created on: 5 Feb 2018
 *      Author: ban115
 */

#ifndef DADASOURCE_H_
#define DADASOURCE_H_

#include "DataSource.h"
#include "CpuTimer.h"
#include "dada_hdu.h"
#include "dada_def.h"
#include "ipcio.h"

class DadaSource: public DataSource {
public:
	DadaSource(int nt, const char* keyname, bool lock);
	virtual ~DadaSource();

	int get_header_int(const char* name);
	double get_header_double(const char* name);
	int get_header_string(const char* name, char* out);


	int npols() {
		return m_npols;
	}
	int nbeams() {
		return m_nbeams;
	}
	int nchans() {
		return m_nchans;
	}
	int nbits() {
		return m_nbits;
	}
	int nants() {
		return 1;
	}
	double tsamp() {
		return m_tsamp;
	}
	double tstart() {
		return m_tstart;
	}
	double fch1() {
		return m_fch1;
	}
	double foff() {
		return m_foff;
	}
	DataOrder data_order() {
		return m_out_data_order;
	}
	size_t current_sample()
	{
	  return m_current_sample;
	}
	size_t read_samples(void** output);
	size_t seek_sample(size_t t);

	void release_buffer();

	// I need to put these here otherwise they get changed
	// At the end of of the constructor, and I have no idea why.
	char* name();

	int dada_key() {
		return m_key;
	}

private:
	void* get_next_buffer(size_t& nt);

	dada_hdu_t* m_hdu;
	char* m_hdr;
    int m_npols;
    int m_nbeams;
    int m_nchans;
    int m_nbits;
    int m_nt;
    double m_tsamp;
    double m_fch1;
    double m_foff;
    double m_tstart;
    bool m_got_buffer;
    DataOrder m_in_data_order;
    DataOrder m_out_data_order;
    uint64_t m_bytes_per_block;
    uint64_t m_blkid;
    size_t m_current_sample;
    void* m_reorder_buffer;
    CpuTimer m_transpose_timer;
    int m_key;

};

#endif /* DADASOURCE_H_ */
