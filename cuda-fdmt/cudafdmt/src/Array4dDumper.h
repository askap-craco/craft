/*
 * Array4dDumper.h
 *
 *  Created on: 15 Apr 2019
 *      Author: ban115
 */

#ifndef ARRAY4DDUMPER_H_
#define ARRAY4DDUMPER_H_

#include "array.h"
#include "DataSource.h"
#include "FreddaParams.h"
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include <string>
using std::string;

class Array4dDumper {

	FILE* m_fout;
	std::string m_name;
	array4d_t& m_target;
	size_t m_num_elements;
	void _fwrite(const void* ptr, size_t size, size_t count);
	bool m_auto_copy;

public:
	Array4dDumper(array4d_t& target, const char* name, FreddaParams& params, int tsamp_mult, bool auto_copy=true);
	virtual ~Array4dDumper();
	void dump();

};


#endif /* ARRAY4DDUMPER_H_ */
