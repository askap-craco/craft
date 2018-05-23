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
#include "InvalidSourceFormat.h"
#include <fcntl.h> // posix_fadvise

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
		throw InvalidSourceFormat();
	}

	// find end of header
	size_t size = fread(m_hdr, sizeof(char), MAX_HDR_SIZE, m_file);
	char* hdr_end = mystrnstr(m_hdr, "HEADER_END", MAX_HDR_SIZE);
	if (hdr_end == NULL) {
		printf("SigprocFile: File %s does not contain HEADER_END\n", filename);
		throw InvalidSourceFormat();
	}
	// TODO: Check it starts with HEADER_START
	m_hdr_nbytes = (size_t)((hdr_end - m_hdr) + strlen("HEADER_END"));
	assert(m_hdr_nbytes < MAX_HDR_SIZE);
	m_hdr[m_hdr_nbytes] = 0;
	seek_sample(0);
	printf("File %s has header %d bytes long. i.e. 0x%x\n", filename, m_hdr_nbytes, m_hdr_nbytes);
	m_nbits = header_int("nbits");
	m_nifs = header_int("nifs");
	m_nchans = header_int("nchans");
	m_fch1 = header_double("fch1");
	m_foff = header_double("foff");
	m_tstart = header_double("tstart");
	m_tsamp = header_double("tsamp");

	// tell_linux we'll be reading sequentailly
	m_fd = fileno(m_file);

	if (posix_fadvise(m_fd, 0, 0, POSIX_FADV_SEQUENTIAL) != 0) {
	  perror("Could not set advice\n");
	  exit(EXIT_FAILURE);
	}
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
	if(fseek(m_file, boff, SEEK_SET) < 0) {
		printf("SigprocFile: Could not seek to offset of file %s\n. Error: %s", m_filename, strerror(errno));
		assert(0);
	}

	m_current_sample = t;
	return boff;
}

void SigprocFile::advise_block(off_t bytes_per_block)
{
  int nblocks = 16;
  off_t offset = ftell(m_file);
  if (posix_fadvise(m_fd, offset, bytes_per_block*nblocks, POSIX_FADV_WILLNEED) != 0) {
    perror("Couln't set advice for next block\n");
    exit(EXIT_FAILURE);
  }

  // tell linux we don't need the stuff we've read
  if (posix_fadvise(m_fd, 0, offset, POSIX_FADV_DONTNEED) != 0) {
    perror("Coulnt set advise for previous data\n");
    exit(EXIT_FAILURE);
  }
}

size_t SigprocFile::read_samples_uint8(size_t nt, uint8_t* output)
{
	// RETURNS TBF ordering. WARNING: This will deeply confuse the output if nbeams != 1
	assert(nifs() == 1); // Otherwise users will be confused. TODO: get sources to tell user what data order is
	assert(m_nbits == 8);
	size_t nreq = nt*m_nifs*m_nchans;
	size_t nelements = fread(output, sizeof(uint8_t), nreq, m_file);
	size_t ont = nelements/m_nifs/m_nchans;
	m_samples_read += ont;
	m_current_sample += nt;
	advise_block(sizeof(uint8_t)*nreq);
	return ont;
}


size_t SigprocFile::read_samples(void** output)
{
	// Usually we use the buffer in the FileSet, so we'll just ignore this for now.
	assert(0);
	return 0;
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
