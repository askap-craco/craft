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

DadaSource::DadaSource(int key, bool lock) {
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
	if (ascii_header_get (header, "HDR_SIZE", "%"PRIu64, &hdr_size) != 1) {
		hdr_size = header_size;
		if (ascii_header_set (header, "HDR_SIZE", "%"PRIu64, hdr_size) < 0) {
			printf("Couldnt set header size");
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
	m_tsamp = get_header_double("TSAMP");
	m_fch1 = get_header_double("FREQ");
	m_foff = get_header_double("BW");
	m_tstart = get_header_double("MJD_START");

}

DadaSource::~DadaSource() {
	if (m_hdr) {
		free(m_hdr);
	}
	if (m_hdu) {
		dada_hdu_destroy(m_hdu);
	}
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
size_t DadaSource::read_samples_uint8(size_t nt, uint8_t* output)
{
	return size_t(0);
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
	return "Hello";
}
