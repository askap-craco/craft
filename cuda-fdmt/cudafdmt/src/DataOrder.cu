/*
 * DataOrder.cpp
 *
 *  Created on: 16 Apr 2018
 *      Author: ban115
 */

#include <stdio.h>
#include <assert.h>
#include "DataOrder.h"


DataOrder data_order_from_string(const char* str)
{
	DataOrder d;
	if (strcmp(str, "TFBP") == 0) {
		d = DataOrder::TFBP;
	} else if (strcmp(str, "FTBP") == 0) {
		d = DataOrder::FTBP;
	} else if (strcmp(str, "BPFT") == 0) {
		d = DataOrder::BPFT;
	} else if (strcmp(str, "BPTF") == 0)  {
		d = DataOrder::BPTF;
	} else {
		printf("Unknown data order: %s\n", str);
		assert(1==0);
		exit(EXIT_FAILURE);
	}

	return d;
}
