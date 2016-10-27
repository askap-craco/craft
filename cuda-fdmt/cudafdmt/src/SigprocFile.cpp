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
	m_filename = new char[strlen(filename)];
	m_filename = strcpy(m_filename, filename);

	if (! m_file) {
		printf("SigprocFile: could not open file: %s - \n",filename, strerror(errno));
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
	size_t hdr_nbytes = (hdr_end - m_hdr) + strlen("HEADER_END");
	assert(hdr_nbytes < MAX_HDR_SIZE);
	m_hdr[hdr_nbytes] = 0;
	fseek(m_file, hdr_nbytes, SEEK_SET);
	printf("File %s has header %d bytes long\n", filename, hdr_nbytes);
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


