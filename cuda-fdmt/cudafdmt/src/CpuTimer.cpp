/*
 * CpuTimer.cpp
 *
 *  Created on: 28 Oct 2016
 *      Author: ban115
 */

#include "CpuTimer.h"
#include <time.h>

timespec diff(timespec start, timespec end)
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

CpuTimer::CpuTimer() {
	// TODO Auto-generated constructor stub

}

CpuTimer::~CpuTimer() {
	// TODO Auto-generated destructor stub
}

//void CpuTimer::start() {
////	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &m_thread_start);
////	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &m_cpu_start);
////	clock_gettime(CLOCK_MONOTONIC, &m_mono_start);
//}

//void CpuTimer::start() {
////	clock_gettime(CLOCK_THREAD_CPUTIME_ID, &m_thread_stop);
////	clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &m_cpu_stop);
////	clock_gettime(CLOCK_MONOTONIC, &m_mono_stop);
//}

