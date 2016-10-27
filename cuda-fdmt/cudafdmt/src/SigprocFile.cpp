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
	char* hdr_end = strstr(m_hdr, "HEADER_END");
	if (hdr_end == NULL) {
		printf("SigprocFile: File %s does not contain HEADER_END\n", filename);
		assert(hdr_end != NULL);
		exit(EXIT_FAILURE);
	}
	// TODO: Check it starts with HEADER_START
	size_t hdr_nbytes = (hdr_end - m_hdr) + strlen("HEADER_END");
	m_hdr[hdr_nbytes] = 0;
	fseek(m_file, hdr_nbytes, SEEK_SET);
	printf("File % has header %d bytes long\n", filename, hdr_nbytes);
}

SigprocFile::~SigprocFile() {
	if(m_file) {
		fclose(m_file);
		delete[] m_filename;
	}
}

float SigprocFile::header(const char* hname) {
	char* hstart = strstr(m_hdr, hname);
	if (hstart == NULL) {
		printf("Could not find header %s in file %s\n", hname, m_filename);
		exit(EXIT_FAILURE);
	}
	char* dstart = (hstart + strlen(hname));
	int* dlen_ptr = (int*) dstart;
	int dlen = (int)*dlen_ptr;
	char* vstart = dstart + sizeof(int);
	float* value_ptr = (float*) vstart;
	float value = (float) *vstart;

	return value;

}


