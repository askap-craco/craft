/*
 * CpuTimer.h
 *
 *  Created on: 28 Oct 2016
 *      Author: ban115
 */

#ifndef CPUTIMER_H_
#define CPUTIMER_H_
#include <time.h>
#include <iostream>

using namespace std;

class CpuTimer {
public:
	CpuTimer();
	virtual ~CpuTimer();

	void start();
	void stop();
	timespec m_thread_start;
	timespec m_proc_start;
	timespec m_mono_start;

	timespec m_thread_stop;
	timespec m_proc_stop;
	timespec m_mono_stop;

	friend ostream &operator<<(ostream & output,  CpuTimer &t);
};

#endif /* CPUTIMER_H_ */
