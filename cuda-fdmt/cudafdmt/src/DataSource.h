/*
 * DataSource.h
 *
 *  Created on: 1 Nov 2016
 *      Author: ban115
 */

#ifndef DATASOURCE_H_
#define DATASOURCE_H_
#include <stdio.h>
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
	virtual size_t seek_sample(size_t t) = 0;
	virtual size_t samples_read() =0;
	virtual char* name() = 0;
	float dm_of_idt(int idt) {

			float nu1 = fch1()/1e3;
			float nu2 = (fch1() + foff()*nchans())/1e3;
			float dm = fabs(idt*tsamp() / 4.15e-3 / (1.0/(nu1*nu1) - 1.0/(nu2*nu2)));

			return dm;

	};


	CpuTimer m_read_timer;



};

#endif /* DATASOURCE_H_ */
