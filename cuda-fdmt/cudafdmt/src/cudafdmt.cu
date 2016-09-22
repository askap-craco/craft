/*
 ============================================================================
 Name        : cudafdmt.cu
 Author      : Keith Bannister <keith.bannister@csiro.au>
 Version     :
 Copyright   : (C) CSIRO 2016
 Description : Compute sum of reciprocals using STL on CPU and Thrust on GPU
 ============================================================================
 */

#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>

#include <thrust/reduce.h>
#include <thrust/device_vector.h>

template <typename T> __host__ __device__  T reciprocal(const T &x)
{
	return ((T)1)/x;
}

template <typename T> class ReciprocalFunctor {
	public:
	__host__ __device__ T operator()(const T &x) {
		return reciprocal(x);
	}
};

template <typename T, class OpClass> T transformAndSumCPU(std::vector<T> data, OpClass op)
{
	std::vector<T> temp(data.size());
	std::transform(data.begin(), data.end(), temp.begin(), op);
	return std::accumulate(temp.begin(), temp.end(), (T)0);
}

template <typename T, class OpClass> T transformAndSumGPU(std::vector<T> data, OpClass op)
{
	thrust::device_vector<T> temp(data.begin(), data.end());
	thrust::transform(temp.begin(), temp.end(), temp.begin(), op);
	return thrust::reduce(temp.begin(), temp.end());
}

template<typename T> void initialize(std::vector<T> &data, unsigned workSize)
{
	/* Initialize the vector */
	for (unsigned i = 0; i < workSize; i++)
		data.push_back( ((T)0.1)*(i+1) );
}

template<typename T> void doCompute(unsigned workSize)
{
	std::vector<T> hostData;

	initialize(hostData, workSize);
	T cpuResults = transformAndSumCPU(hostData, ReciprocalFunctor<T>());
	T gpuResults = transformAndSumGPU(hostData, ReciprocalFunctor<T>());
	std::cout<<"transformAndSumCPU = "<<cpuResults<<std::endl;
	std::cout<<"transformAndSumGPU = "<<gpuResults<<std::endl;
}

