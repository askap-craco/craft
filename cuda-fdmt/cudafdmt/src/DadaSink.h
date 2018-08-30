/*
 * DadaSink.h
 *
 *  Created on: 29 Aug 2018
 *      Author: ban115
 */

#ifndef DADASINK_H_
#define DADASINK_H_

#include "DataSource.h"
#include "dada_hdu.h"
#include "dada_def.h"
#include "ipcio.h"

class DadaSink {

	dada_hdu_t* m_hdu;
	uint64_t m_blkid;
	uint64_t m_data_block_size;
	void* m_current_block;

public:
	DadaSink(DataSource& source, int key, char* hdr,
			int npol_out, int nbeams_out, int nt);
	virtual ~DadaSink();
	inline void* current_block() { return m_current_block; }
	void* open_block();
	void close_block();

};

#endif /* DADASINK_H_ */

