/*
 * CudaTimer.h
 *
 *  Created on: 4 Oct 2016
 *      Author: ban115
 */

#ifndef CUDATIMER_H_
#define CUDATIMER_H_

#include <iostream>
#include <cuda.h>
#include <cuda_runtime.h>

using namespace std;

class CudaTimer {
public:
	CudaTimer(cudaStream_t stream = 0);
	virtual ~CudaTimer();

	void start();
	void stop();
	void sync_start();
	void sync_stop();
	float get_elapsed_time();

	friend ostream &operator<<(ostream & output,  CudaTimer &t)
	{
		output << "total: " << t.m_total_time
				<< " average: " << t.m_total_time/t.m_ncalls
				<< " last: " << t.get_elapsed_time() << " ms"
				<< " ncalls:" <<	t.m_ncalls;

		return output;
	}

private:
	cudaEvent_t m_start;
	cudaEvent_t m_stop;
	cudaStream_t m_stream;
	float m_total_time;
	float m_ncalls;
};

#endif /* CUDATIMER_H_ */
