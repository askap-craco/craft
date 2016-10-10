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
	void stop(bool sync=true);
	void sync_start();
	void sync_stop();
	float get_elapsed_time();

	friend ostream &operator<<(ostream & output,  CudaTimer &t)
	{
		output << t.get_elapsed_time() << " ms";
		return output;
	}

private:
	cudaEvent_t m_start;
	cudaEvent_t m_stop;
	cudaStream_t m_stream;
};

#endif /* CUDATIMER_H_ */
