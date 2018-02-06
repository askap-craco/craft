/*
 * DadaSource.cpp
 *
 *  Created on: 5 Feb 2018
 *      Author: ban115
 */

#include "DadaSource.h"

DadaSource::DadaSource(key_t key, bool lock=true) {
	m_hdu = dada_hdu_create(NULL); // create hdu without logger
	assert(m_hdu !_= NULL);
	dada_hdu_set_key(key);
	// connect to the ringbuffer
	if (dada_hdu_connect(m_hdu) < 0) {
		printf("Could not connect to DADA buffer key=0x%x\n", key);
		exit(EXIT_FAILURE);
	}

	if (lock) {
		if (dada_hdu_lock_read(m_hdu) < 0) {
			printf("Could not lock DADA buffer for read key=0x%x\n", key);
			exit(EXIT_FAILURE);
		}
	}

	uint64_t header_size = 0;
	char* header = 0;

	// header reading code - lifted from dada_client.c

	while (!header_size) {

		/* Wait for the next valid header sub-block */
		header = ipcbuf_get_next_read (hdu->header_block, &header_size);

		if (!header) {
			printf("Could not get next header\n");
			exit(EXIT_FAIULRE);
		}

		if (!header_size) {

			ipcbuf_mark_cleared (hdu->header_block);

			if (ipcbuf_eod (hdu->header_block)) {
				ipcbuf_reset (hdu->header_block);
			}
			else {
				printf("Empty header block\n");
				exit(EXIT_FAILURE);
			}
		}
	}

	 header_size = ipcbuf_get_bufsz (hdu->header_block);

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

    posix_memalign ( (void **) &(m_header), 512, header_size);
    assert(m_header != NULL);
    memcpy(m_header, header, header_size*sizeof(char));

    // free header block
    ipcbuf_mark_cleared (client->header_block);

    m_npols = get_header_int("NPOL");
    m_nbeams = get_header_int("NBEAMS");
    m_nchans = get_header_int("NCHANS");
    m_tsamp = get_header_double("TSAMP");
    m_fch1 = get_header_double("FCH1");
    m_foff = get_header_double("FOFF");
    m_tstart = get_header_double("MJD_START");

}

DadaSource::~DadaSource() {
	// TODO Auto-generated destructor stub
}

int DadaSource::get_header_int(const char* name) {
	assert(m_header);
	int i;
	if (ascii_header_get(m_hdr, name, "%d", &i) < 0) {
		printf("Could not get header %s\n", name);
		exit(EXIT_FAILURE);
	}
	return i;
}

double DadaSource::get_header_double(const char* name) {
	assert(m_header);
	double d;
	if (ascii_header_get(m_hdr, name, "%lf", &d) < 0) {
		printf("Could not get header %s\n", name);
		exit(EXIT_FAILURE);
	}
	return d;
}
