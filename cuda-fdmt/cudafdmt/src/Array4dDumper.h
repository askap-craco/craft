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


class Array4dDumper {

	FILE* m_fout;

public:
	Array4dDumper(array4d_t& target, const char* name, FreddaParams& params);
	virtual ~Array4dDumper();
};


#endif /* ARRAY4DDUMPER_H_ */
