/*
 * DadaSource.h
 *
 *  Created on: 5 Feb 2018
 *      Author: ban115
 */

#ifndef DADASOURCE_H_
#define DADASOURCE_H_

#include "DataSource.h"
#include "dada_hdu.h"
#include "dada_def.h"
#include "ipcio.h"

class DadaSource: public fdmt::DataSource {
public:
	DadaSource(key_t key, bool lock);
	virtual ~DadaSource();

	int get_header_int(const char* name);
	double get_header_double(const char* name);

	int npols() {
		return m_npols;
	}
	int nbeams() {
		return m_nbeams;
	}
	int nchans() {
		return m_nchans;
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

private:
	dada_hdu_t* m_hdu;
	char* m_hdr;
    int m_npols;
    int m_nbeams;
    int m_nchans;
    double m_tsamp;
    double m_fch1;
    double m_foff;
    double m_tstart;

};

#endif /* DADASOURCE_H_ */
