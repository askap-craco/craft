/*
 * DataSource.h
 *
 *  Created on: 1 Nov 2016
 *      Author: ban115
 */

#ifndef DATASOURCE_H_
#define DATASOURCE_H_
#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include "CpuTimer.h"
#include "DataOrder.h"


class DataSource {
public:
	DataSource();
	virtual ~DataSource();

	virtual double fch1() = 0;
	virtual double foff() = 0;
	virtual double tstart() = 0;
	virtual double tsamp() = 0;
	virtual int nants() = 0;
	virtual int nbeams() = 0;
	virtual int npols() = 0;
	virtual int nchans() = 0;
	virtual int nbits() = 0;
	virtual DataOrder data_order() = 0;

	virtual size_t read_samples(void** output) = 0;
	virtual size_t read_samples_ant(void** output, int iant) {
		return read_samples(output);
	}

	virtual size_t seek_sample(size_t t) = 0;
	virtual size_t current_sample() = 0;
	virtual double current_mjd();
	
	virtual char* name() = 0;
	virtual const char* antenna_name() { return ""; };
	float dm_of_idt(int idt) {

			float nu1 = fch1()/1e3;
			float nu2 = (fch1() + foff()*nchans())/1e3;
			float dm = fabs(idt*tsamp() / 4.15e-3 / (1.0/(nu1*nu1) - 1.0/(nu2*nu2)));

			return dm;

	};

	int64_t current_sample_relative_to(double mjd) {
		double mjddiff = mjd - current_mjd();
		int64_t sampdiff = (int64_t)round(mjddiff*86400./tsamp());
		return sampdiff;
	}


	CpuTimer m_read_timer;



};

#endif /* DATASOURCE_H_ */
