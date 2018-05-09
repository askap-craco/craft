/*
 * DadaSource.cpp
 *
 *  Created on: 5 Feb 2018
 *      Author: ban115
 */

#include "DadaSource.h"
#include "ascii_header.h"
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "InvalidSourceFormat.h"

DadaSource::DadaSource(int nt, int key, bool lock) {
	m_got_buffer = false;
	m_hdu = dada_hdu_create(NULL); // create hdu without logger
	assert(m_hdu != NULL);
	dada_hdu_set_key(m_hdu, (key_t)key);
	// connect to the ringbuffer
	if (dada_hdu_connect(m_hdu) < 0) {
		printf("Could not connect to DADA buffer key=0x%x\n", key);
		throw InvalidSourceFormat();
	}

	if (lock) {
		if (dada_hdu_lock_read(m_hdu) < 0) {
			printf("Could not lock DADA buffer for read key=0x%x\n", key);
			exit(EXIT_FAILURE);
		}
	}

	uint64_t header_size = 0, hdr_size = 0;
	char* header = 0;

	// header reading code - lifted from dada_client.c
	while (!header_size) {

		/* Wait for the next valid header sub-block */
		header = ipcbuf_get_next_read (m_hdu->header_block, &header_size);

		if (!header) {
			printf("Could not get next header\n");
			exit(EXIT_FAILURE);
		}

		if (!header_size) {
			ipcbuf_mark_cleared (m_hdu->header_block);

			if (ipcbuf_eod (m_hdu->header_block)) {
				ipcbuf_reset (m_hdu->header_block);
			}
			else {
				printf("Empty header block\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	header_size = ipcbuf_get_bufsz (m_hdu->header_block);

	/* Check that header is of advertised size */
	if (ascii_header_get (header, "HDR_SIZE", "%" PRIu64, &hdr_size) != 1) {
		hdr_size = header_size;
		if (ascii_header_set (header, "HDR_SIZE", "%" PRIu64, hdr_size) < 0) {
			printf("Couldn't set header size");
			exit(EXIT_FAILURE);
		}
	}

	if (hdr_size < header_size)
		header_size = hdr_size;

	posix_memalign ( (void **) &(m_hdr), 512, header_size);
	assert(m_hdr != NULL);
	memcpy(m_hdr, header, header_size*sizeof(char));

	// free header block
	ipcbuf_mark_cleared (m_hdu->header_block);

	m_npols = get_header_int("NPOL");
	m_nbeams = get_header_int("NBEAM");
	m_nchans = get_header_int("NCHAN");
	m_nbits = get_header_int("NBIT");
	m_tsamp = get_header_double("TSAMP");
	m_fch1 = get_header_double("FREQ");
	m_foff = get_header_double("BW");
	m_tstart = get_header_double("MJD_START");
	m_bytes_per_block = npols()*nbeams()*nchans()*nbits()*nt/8;

	char order_str[256];
	get_header_string("DORDER",order_str);
	m_reorder_buffer = NULL;
	m_in_data_order = data_order_from_string(order_str);
	if (m_in_data_order == DataOrder::TFBP) {
		m_reorder_buffer = malloc(m_bytes_per_block);
		assert(m_reorder_buffer);
	} else {
		printf("Unknown order: %s\b", order_str);
		exit(EXIT_FAILURE);
	}

	// output data order
	//m_out_data_order = DataOrder::BPTF;
	m_out_data_order = DataOrder::TFBP;
	assert(m_in_data_order == DataOrder::TFBP);
}

DadaSource::~DadaSource() {
	if (m_hdr) {
		free(m_hdr);
	}
	if (m_hdu) {
		dada_hdu_destroy(m_hdu);
	}
	if (m_reorder_buffer) {
		free(m_reorder_buffer);
	}
	cout << "DADA Transpose timing " << m_transpose_timer << endl;
}

int DadaSource::get_header_int(const char* name) {
	assert(m_hdr);
	int i;
	if (ascii_header_get(m_hdr, name, "%d", &i) < 0) {
		printf("Could not get header %s\n", name);
		exit(EXIT_FAILURE);
	}
	return i;
}

double DadaSource::get_header_double(const char* name) {
	assert(m_hdr);
	double d;
	if (ascii_header_get(m_hdr, name, "%lf", &d) < 0) {
		printf("Could not get header %s\n", name);
		exit(EXIT_FAILURE);
	}
	return d;
}

int DadaSource::get_header_string(const char* name, char* out) {
	assert(m_hdr);
	if (ascii_header_get(m_hdr, name, "%s", out) < 0) {
		printf("Could not get header %s\n",name);
		exit(EXIT_FAILURE);
	}

	return 0;
}

size_t DadaSource::read_samples(void** output)
{
	m_read_timer.start();
	// TODO: get main loop to release as soon as it's finished, rather than on next read, but it's too too much typing.
	if (m_got_buffer) {
		release_buffer();
	}
	uint64_t nbytes;
	char* ptr = ipcio_open_block_read(m_hdu->data_block, &nbytes, &m_blkid);
	m_read_timer.stop();

	m_got_buffer = true;
	// TODO check expected nbytes
	//m_bytes_per_block = npols()*nbeams()*nchans()*nbits()*nt/8;
	size_t nt = nbytes/(npols()*nbeams()*nchans()*nbits()/8);
	assert(m_in_data_order == DataOrder::TFBP);


	m_transpose_timer.start();
	if (m_in_data_order == m_out_data_order) {
		*output = (void*)ptr;
	} else if (m_in_data_order == DataOrder::TFBP && m_out_data_order == DataOrder::BPTF) {
		*output = m_reorder_buffer;
		assert(nbits() == 32);
		float *inp, *outp;
		inp = (float*) ptr;
		outp = (float*) *output;
		for(int t = 0; t < nt; ++t) {
			for (int f = 0; f < nchans(); ++f) {
				for (int b = 0; b < nbeams(); ++b) {
					for(int p = 0; p < npols(); ++p) {
						int outidx = f + nchans()*(t + nt*(p + npols()*b));
						assert(outidx >= 0);
						assert(outidx < nchans()*npols()*nbeams()*nt);
						outp[outidx] = *inp;
						++inp;
					}
				}
			}
		}
	} else {
		printf("Invalid ordering\n");
		assert(1==0);
		exit(EXIT_FAILURE);
	}
	m_transpose_timer.stop();
	return nt;
}

void DadaSource::release_buffer() {
	if(m_got_buffer) {
		if (ipcio_close_block_read(m_hdu->data_block, m_bytes_per_block)) {
			printf("Could not mark cleared");
			exit(1);
		}

		m_got_buffer = false;
	}
}
size_t DadaSource::seek_sample(size_t t)
{
	return size_t(0);
}
size_t DadaSource::samples_read()
{
	return size_t(0);
}
char* DadaSource::name()
{
	return (char*)"Hello";
}
