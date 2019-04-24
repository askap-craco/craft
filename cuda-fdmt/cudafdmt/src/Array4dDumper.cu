/*
 * Array4dDumper - writes an array4d to disk in DADA format
 *
 *  Created on: 15 April 2019
 *      Author: ban115
 */


#include "Array4dDumper.h"
#include "ascii_header.h"


using std::string;
const int HDR_SIZE = 16384;

Array4dDumper::Array4dDumper(array4d_t& target, const char* name, FreddaParams& params, bool auto_copy)
: m_target(target), m_name(name), m_auto_copy(auto_copy) {
	std::string fname(name);
	fname += ".dada";
	m_num_elements = target.nw*target.nx*target.ny*target.nz;
	assert(m_num_elements > 0);
	m_fout = fopen(fname.c_str(), "w+");
	if (m_fout == NULL) {
		perror("Error opening dumpfile:");
		exit(EXIT_FAILURE);
	}

	char header_buf[HDR_SIZE];
	bzero(header_buf, HDR_SIZE);
	ascii_header_set(header_buf, "HDR_VERSION", "%s", "1.0");
	ascii_header_set(header_buf, "HDR_SIZE", "%d", HDR_SIZE);
	ascii_header_set(header_buf, "INSTRUMENT", "%s", "FREDDA");
	ascii_header_set(header_buf, "FREDDA_VERSION", "%s", VERSION);
	ascii_header_set(header_buf, "DATA_NAME", "%s", name);
	ascii_header_set(header_buf, "NW", "%d", target.nw);
	ascii_header_set(header_buf, "NX", "%d", target.nx);
	ascii_header_set(header_buf, "NY", "%d", target.ny);
	ascii_header_set(header_buf, "NZ", "%d", target.nz);
	ascii_header_set(header_buf, "NUM_ELEMENTS", "%d", m_num_elements);
	ascii_header_set(header_buf, "FILE_NAME", "%s", fname.c_str());
	ascii_header_set(header_buf, "SHAPE", "%d,%d,%d,%d", target.nw,target.nx,target.ny,target.nz);
	ascii_header_set(header_buf, "DTYPE", "%s", "<f4"); // 32 bit little endian floats
	ascii_header_set(header_buf, "TSAMP", "%0.12f", params.source->tsamp()*params.nt*params.num_rescale_blocks); // tsamp unchanged from source
	params.to_dada(header_buf);
	_fwrite(header_buf, sizeof(char), HDR_SIZE);
	fflush(m_fout);
}

Array4dDumper::~Array4dDumper() {
	assert(m_fout != NULL);
	fflush(m_fout);
	fclose(m_fout);
	m_fout = NULL;
}

void Array4dDumper::_fwrite(const void* ptr, size_t size, size_t count) {
	size_t nelements = fwrite(ptr, size, count, m_fout);
	if (nelements != count) {
		perror("Error writing to dumpfile");
		exit(EXIT_FAILURE);
		// Now what do we do? Keep going? Or fail? or close?
	}
}

void Array4dDumper::dump() {
	if (m_auto_copy) {
		array4d_copy_to_host(&m_target);
	}
	int nzero = 0;
	for(int i = 0; i < m_num_elements; ++i) {
		if (m_target.d[i] == 0.0f) {
			nzero += 1;
		}
	}
	printf("%s contained %d zeros\n", m_name.c_str(), nzero);
	_fwrite(m_target.d, sizeof(float), m_num_elements);
}
