/*
 * CpuTimer.cpp
 *
 *  Created on: 28 Oct 2016
 *      Author: ban115
 */

#include "CpuTimer.h"
#include <time.h>
#include <sys/time.h>
#include <stdio.h>
#include <iostream>

#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
//#define CLOCK_THREAD_CPUTIME_ID CALENDAR_CLOCK
//#define CLOCK_PROCESS_CPUTIME_ID CALENDAR_CLOCK
//#define CLOCK_MONOTONIC CALENDAR_CLOCK
//typedef long clockid_t;
#endif

timespec diff(timespec& start, timespec& end)
{
	timespec temp;
	if ((end.tv_nsec-start.tv_nsec)<0) {
		temp.tv_sec = end.tv_sec-start.tv_sec-1;
		temp.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
	} else {
		temp.tv_sec = end.tv_sec-start.tv_sec;
		temp.tv_nsec = end.tv_nsec-start.tv_nsec;
	}
	return temp;
}

double ts2secs(timespec& t) {
	double v = (double)t.tv_sec;
	v += ((double)t.tv_nsec)/1e9;

	return v;
}

// From here:https://gist.github.com/jbenet/1087739
void current_time(clockid_t clk_id, struct timespec *ts) {

#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
  clock_serv_t cclock;
  mach_timespec_t mts;
  host_get_clock_service(mach_host_self(), clk_id, &cclock);
  clock_get_time(cclock, &mts);
  mach_port_deallocate(mach_task_self(), cclock);
  ts->tv_sec = mts.tv_sec;
  ts->tv_nsec = mts.tv_nsec;
#else
  clock_gettime(clk_id, ts);
#endif

}

CpuTimer::CpuTimer() {
	m_ncalls = 0;
	m_cpu_total = 0;
	m_thread_total = 0;
	m_mono_total = 0;
}

CpuTimer::~CpuTimer() {
	// TODO Auto-generated destructor stub
}

void CpuTimer::start() {
     	current_time(CLOCK_THREAD_CPUTIME_ID, &m_thread_start);
	current_time(CLOCK_PROCESS_CPUTIME_ID, &m_cpu_start);
	current_time(CLOCK_MONOTONIC, &m_mono_start);
}

void CpuTimer::stop() {
	current_time(CLOCK_THREAD_CPUTIME_ID, &m_thread_stop);
	current_time(CLOCK_PROCESS_CPUTIME_ID, &m_cpu_stop);
	current_time(CLOCK_MONOTONIC, &m_mono_stop);

	m_cpu_total += ts2secs(m_cpu_stop) - ts2secs(m_cpu_start);
	m_thread_total += ts2secs(m_thread_stop) - ts2secs(m_thread_start);
	m_mono_total += ts2secs(m_mono_stop) - ts2secs(m_mono_start);
	m_ncalls += 1;

}


