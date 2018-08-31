/*
 * DadaSink.cpp
 *
 *  Created on: 29 Aug 2018
 *      Author: ban115
 */

#include "DadaSink.h"
#include "rescale.h"
#include "ascii_header.h"


DadaSink::DadaSink(DataSource& source, int key, char* hdr,
		int npol_out, int nbeams_out, int nt
) {
	m_hdu = dada_hdu_create(NULL);
	dada_hdu_set_key(m_hdu, key);
	if (dada_hdu_connect(m_hdu) < 0) {
		printf("Could not connect to output DADA key 0x%x", key);
		exit(EXIT_FAILURE);
	}
	if (dada_hdu_lock_write(m_hdu) < 0) {
		printf("Could not lock write DADA key 0x%x", key);
		exit(EXIT_FAILURE);
	}

	uint64_t header_size = ipcbuf_get_bufsz(m_hdu->header_block);
	m_data_block_size = ipcbuf_get_bufsz(&m_hdu->data_block->buf);
	char* header_buf = ipcbuf_get_next_write(m_hdu->header_block);
	if (hdr != NULL) { // copy header to output
		int hdrlen = strlen(hdr);
		assert(hdrlen < header_size);
		memcpy(header_buf, hdr, hdrlen*sizeof(char));
	}
	// rescaling always puts the minimum frequency as the channel 0
	double out_freq;
	if (source.foff() > 0) {
		out_freq = source.fch1();
	} else {
		out_freq = source.fch1() + source.nchans()*source.foff();
	}
	assert(out_freq <= source.fch1());
	double out_foff = abs(source.foff());

	uint64_t expected_block_size = npol_out*nbeams_out*source.nchans()*nt*sizeof(rescale_dtype);
	if (m_data_block_size != expected_block_size) {
		printf("DADA output is block size is %d but expected %d\n",
				m_data_block_size, expected_block_size);
		exit(EXIT_FAILURE);
	}
	// set some header parameters - probably could set more, but it's a pain. Sheesh metadata is trickly
	ascii_header_set(header_buf, "HDR_VERSION", "%s", "1.0");
	ascii_header_set(header_buf, "INSTRUMENT", "%s", "FREDDA_ICS");
	ascii_header_set(header_buf, "FREDDA_VERSION", "%s", VERSION);
	ascii_header_set(header_buf, "NPOL", "%d", npol_out); // polarisations summed
	ascii_header_set(header_buf, "NBEAMS", "%d", nbeams_out); // polarisations summed
	ascii_header_set(header_buf, "NCHAN", "%d", source.nchans()); // polarisations summed
	ascii_header_set(header_buf, "NBIT", "%d", sizeof(rescale_dtype)*8); // floating point numbers
	ascii_header_set(header_buf, "NT", "%d", nt);
	ascii_header_set(header_buf, "TSAMP", "%0.12f", source.tsamp()); // tsamp unchanged from source
	ascii_header_set(header_buf, "FREQ", "%0.12f", out_freq); // first channel
	ascii_header_set(header_buf, "BW", "%0.12f", out_foff); // Frequency offset - different from source as rescaler sets first frequency to bottom always
	ascii_header_set(header_buf, "MJD_START", "%0.12f", source.current_mjd()); // mjd of first sample
	ascii_header_set(header_buf, "DORDER", "%s", "BPFT");
	ascii_header_set(header_buf, "HDR_SIZE", "%d", header_size);
	ascii_header_set(header_buf, "DATA_TYPE", "%s", "CRAFT");

	if (ipcbuf_mark_filled(m_hdu->header_block, header_size) < 0) {
		printf("Could not mark filled dada sink header block\n");
		exit(EXIT_FAILURE);
	}

	// make memory as pinned so cuda will transfer it more efficiently.
	size_t nbufs = ipcbuf_get_nbufs(&m_hdu->data_block->buf);
	size_t blck_size = ipcbuf_get_bufsz(&m_hdu->data_block->buf);
	for(int b = 0; b < nbufs; b++) {
		gpuErrchk(cudaHostRegister(m_hdu->data_block->buf.buffer[b], blck_size, cudaHostRegisterDefault));
	}
	m_current_block = NULL;
}

DadaSink::~DadaSink() {
	if (m_current_block != NULL) {
		close_block();
	}
	if (m_hdu != NULL) {
		dada_hdu_unlock_write(m_hdu); // Could check for error i.e. < 0, - but whats thie point?
		dada_hdu_disconnect(m_hdu);
	}
}

void* DadaSink::open_block() {
	assert(m_current_block == NULL);
	m_current_block = ipcio_open_block_write(m_hdu->data_block, &m_blkid);
	if (m_current_block == NULL) {
		printf("DADA sink could not open block\n");
		exit(EXIT_FAILURE);
	}
	return m_current_block;
}

void DadaSink::close_block() {
	assert(m_current_block != NULL);
	if (ipcio_close_block_write(m_hdu->data_block, m_data_block_size) < 0) {
		printf("DADA sink could not close block %d\n", m_blkid);
		exit(EXIT_FAILURE);
	}
	m_current_block = NULL;
}
