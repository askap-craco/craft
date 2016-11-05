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
	friend ostream &operator<<(ostream & output,  CpuTimer &t)
	{
		output << "Wall total: " << t.m_mono_total
				<< "s Wall average: " << t.m_mono_total/t.m_ncalls
				<< " CPU total: " << t.m_cpu_total << "s "
				<< " CPU average:" <<	t.m_cpu_total/t.m_ncalls;

		return output;
	}
private:
	timespec m_thread_start;
	timespec m_cpu_start;
	timespec m_mono_start;

	timespec m_thread_stop;
	timespec m_cpu_stop;
	timespec m_mono_stop;
	double m_cpu_total;
	double m_thread_total;
	double m_mono_total;
	double m_ncalls;
};

#endif /* CPUTIMER_H_ */
