/*
 * Array4dDumper - writes an array4d to disk in DADA format
 *
 *  Created on: 15 April 2019
 *      Author: ban115
 */

#include <string>

#include "Array4dDumper.h"
#include "ascii_header.h"


using std::string;
const int HDR_SIZE = 16384;

Array4dDumper::Array4dDumper(array4d_t& target, const char* name, FreddaParams& params) {
	std::string fname(name);
	fname += ".dada";
	m_fout = fopen(fname.c_str(), "w+");
	if (m_fout == NULL) {
		perror("Could not open dumpfile:");
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
	ascii_header_set(header_buf, "SHAPE", "%d,%d,%d,%d", target.nw,target.nx,target.ny,target.nz);
	ascii_header_set(header_buf, "DTYPE", "%s", "<f32"); // 32 bit little endian floats
	ascii_header_set(header_buf, "TSAMP", "%0.12f", params.source->tsamp()*params.nt*params.num_rescale_blocks); // tsamp unchanged from source
	params.to_dada(header_buf);
}

Array4dDumper::~Array4dDumper() {
	assert(m_fout != NULL);
	fclose(m_fout);
	m_fout = NULL;

}




