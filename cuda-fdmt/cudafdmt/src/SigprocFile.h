/*
 * SigprocFile.h
 *
 *  Created on: 27 Oct 2016
 *      Author: ban115
 */

#ifndef SIGPROCFILE_H_
#define SIGPROCFILE_H_

#include <stdio.h>
#include <stdlib.h>

const size_t MAX_HDR_SIZE = 4096;

class SigprocFile {
public:
	SigprocFile(const char* filename);
	virtual ~SigprocFile();
	const char* header_find(const char* hname) const;
	int header_int(const char* hname) const;
	double header_double(const char* hname) const;


private:
	FILE* m_file;
	char* m_filename;
	char m_hdr[MAX_HDR_SIZE];

};

#endif /* SIGPROCFILE_H_ */
