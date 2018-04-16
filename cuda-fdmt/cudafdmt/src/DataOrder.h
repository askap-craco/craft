/*
 * DataOrder.h
 *
 *  Created on: 16 Apr 2018
 *      Author: ban115
 */

#ifndef DATAORDER_H_
#define DATAORDER_H_

enum class DataOrder { TFBP, FTBP, BPFT, BPTF };

DataOrder data_order_from_string(const char* str);

#endif /* DATAORDER_H_ */
