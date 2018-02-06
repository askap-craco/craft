/*
 * InvalidSourceFormat.h
 *
 *  Created on: 6 Feb 2018
 *      Author: ban115
 */

#ifndef INVALIDSOURCEFORMAT_H_
#define INVALIDSOURCEFORMAT_H_

#include <stdexcept>

class InvalidSourceFormat : public std::runtime_error {
public:
	InvalidSourceFormat() : std::runtime_error("InvalidSourceFormat") {};
	InvalidSourceFormat(const std::string& what_arg) : std::runtime_error(what_arg)  {} ;
	InvalidSourceFormat(const char* what_arg) : std::runtime_error(what_arg) {};

};

#endif /* INVALIDSOURCEFORMAT_H_ */
