/*
 * SigprocFile.cpp
 *
 *  Created on: 27 Oct 2016
 *      Author: ban115
 */

#include "SigprocFile.h"

#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <sys/types.h>
#include <stdint.h>
#include <math.h>

/* Same as strstr but goes through *all* the string - even if it contains nulls
 *
 */
char* mystrnstr(const char* s1, const char* s2, size_t n)
{
	size_t s2_len = strlen(s2);
	for (size_t i = 0; i < n - s2_len; i++) {
		if(strncmp(s1+i, s2, s2_len) == 0) {
			return (char*) ( s1 + i);
		}
	}

	return NULL;
}


SigprocFile::SigprocFile(const char* filename) {
	m_file = fopen(filename, "r");
	m_filename = new char[strlen(filename) + 1];
	m_filename = strcpy(m_filename, filename);
	m_samples_read = 0;

	if (! m_file) {
		printf("SigprocFile: could not open file: %s - %s\n",filename, strerror(errno));
		assert(m_file);
		exit(EXIT_FAILURE);
	}

	// find end of header
	size_t size = fread(m_hdr, sizeof(char), MAX_HDR_SIZE, m_file);
	char* hdr_end = mystrnstr(m_hdr, "HEADER_END", MAX_HDR_SIZE);
	if (hdr_end == NULL) {
		printf("SigprocFile: File %s does not contain HEADER_END\n", filename);
		assert(hdr_end != NULL);
		exit(EXIT_FAILURE);
	}
	// TODO: Check it starts with HEADER_START
	m_hdr_nbytes = (size_t)((hdr_end - m_hdr) + strlen("HEADER_END") + 1);
	assert(m_hdr_nbytes < MAX_HDR_SIZE);
	m_hdr[m_hdr_nbytes] = 0;
	seek_sample(0);
	printf("File %s has header %d bytes long\n", filename, m_hdr_nbytes);
	m_nbits = header_int("nbits");
	m_nifs = header_int("nifs");
	m_nchans = header_int("nchans");
	m_fch1 = header_double("fch1");
	m_foff = header_double("foff");
	m_tstart = header_double("tstart");
	m_tsamp = header_double("tsamp");
}

SigprocFile::~SigprocFile() {
	if(m_file) {
		fclose(m_file);
		delete[] m_filename;
	}
}

const char* SigprocFile::header_find(const char* hname) const
{
	char* hstart = mystrnstr(m_hdr, hname, MAX_HDR_SIZE);
	if (hstart == NULL) {
		printf("Could not find header %s in file %s\n", hname, m_filename);
		exit(EXIT_FAILURE);
	}

	int hlen = *((int*)(hstart - sizeof(int)));
	if (hlen != strlen(hname)) {
		printf("Could not find header %s n file %s (but found a substring)\n", hname, m_filename);
		exit(EXIT_FAILURE);
	}

	char* dstart = (hstart + strlen(hname));
	return dstart;
}
double SigprocFile::header_double(const char* hname) const {

	double* vstart = (double*) header_find(hname);
	double value = *vstart;
	return value;
}

int SigprocFile::header_int(const char* hname) const {
	int* vstart = (int*) header_find(hname);
	int value = *vstart;
	return value;
}

size_t SigprocFile::seek_sample(size_t t)
{
	size_t boff = t*nifs()*nchans() + m_hdr_nbytes;
	//printf("Seek t=%d nifs=%d nchans%d m_hdr_nbytes %d\n", t, nifs(), nchans(), m_hdr_nbytes);
	printf("Seek nifs=%d\n", nifs());
	if(fseek(m_file, boff, SEEK_SET) < 0) {
		printf("SigprocFile: Could not seek to offset of file %s\n. Error: %s", m_filename, strerror(errno));
		assert(0);
	}
	return boff;
}

size_t SigprocFile::read_samples_uint8(size_t nt, uint8_t* output)
{
	assert(m_nbits == 8);
	size_t nreq = nt*m_nifs*m_nchans;
	size_t nelements = fread(output, sizeof(uint8_t), nreq, m_file);
	size_t ont = nelements/m_nifs/m_nchans;
	m_samples_read += ont;
	return ont;
}

double SigprocFile::last_sample_elapsed_seconds()
{
	double toff_sec = ((double) m_samples_read) * m_tsamp;
	return toff_sec;
}

double SigprocFile::last_sample_mjd()
{
	double mjd = m_tstart + last_sample_elapsed_seconds()/86400.0;
	return mjd;
}

float SigprocFile::dm_of_idt(int idt)
{
	float nu1 = m_fch1/1e3;
	float nu2 = (m_fch1 + m_foff*m_nchans)/1e3;
	float dm = fabs(idt*m_tsamp / 4.15e-3 / (1.0/(nu1*nu1) - 1.0/(nu2*nu2)));

	return dm;
}

