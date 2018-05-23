/*
 * DataSource.cpp
 *
 *  Created on: 1 Nov 2016
 *      Author: ban115
 */

#include "DataSource.h"


DataSource::DataSource() {
	// TODO Auto-generated constructor stub

}

DataSource::~DataSource() {
	// TODO Auto-generated destructor stub
}

double DataSource::current_mjd() {
  double mjd = tstart() + current_sample()*tsamp()/86400.0;
  
  return mjd;
}
