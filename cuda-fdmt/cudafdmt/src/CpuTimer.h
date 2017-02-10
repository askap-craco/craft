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
		output << "Wall total: " << t.wall_total()
				<< "s Wall average: " << t.wall_average()
				<< " CPU total: " << t.cpu_total() << "s "
				<< " CPU average:" <<	t.cpu_average();

		return output;
	}

	double wall_total() {
		return m_mono_total;
	}

	double wall_average() {
		return wall_total()/m_ncalls;
	}
	double cpu_total() {
		return m_cpu_total;
	}
	double cpu_average() {
		return cpu_total()/m_ncalls;
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
