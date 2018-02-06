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

namespace fdmt {

class DataSource {
public:
	DataSource();
	virtual ~DataSource();

	virtual double fch1();
	virtual double foff();
	virtual double tstart();
	virtual double tsamp();
	virtual int nants();
	virtual int nbeams();
	virtual int npols();
	virtual int nchans();

	virtual size_t read_samples_uint8(size_t nt, uint8_t* output) = 0;
	virtual size_t seek_sample(size_t t) = 0;
	virtual size_t samples_read() =0;
	virtual char* name() = 0;
	virtual float dm_of_idt(int idt);


};

} /* namespace fdmt */

#endif /* DATASOURCE_H_ */
