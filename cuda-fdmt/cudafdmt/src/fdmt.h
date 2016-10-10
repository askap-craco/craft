/**
 * FDMT - based on the
 *
 * Based on FDMT.py
 * Copyright (c) 2014, Barak Zackay (Weizmann Institute of Science)
 * All rights reserved.
 * http://arXiv.org/abs/1411.5373
 */

#ifndef _FDMT_H
#define _FDMT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <vector>
#include "array.h"
#include "array2d.h"

//#define DUMP_STATE 1

#define MAX_ITER 16

using namespace std;

class FdmtIteration
{
public:
	FdmtIteration::FdmtIteration()
	{

	}

	FdmtIteration::FdmtIteration(int nbeam, int nf, int ndt, int nt) {
		//dt_middle_larger = vector<Array1d<int> >();
		//dt_middle_index = vector<Array1d<int> >();
		//dt_rest_index = vector<Array1d<int> >();
		state_shape.w = nbeam;
		state_shape.x = nf;
		state_shape.y = ndt;
		state_shape.z = nt;
	}

	vector<Array1d<int> > dt_middle_larger;
	vector<Array1d<int> > dt_middle_index;
	vector<Array1d<int> > dt_rest_index;
	vector<int> delta_ts;
	coord4_t state_shape;
	void add_subband(int delta_t) {
		dt_middle_larger.push_back(Array1d<int>(delta_t));
		dt_middle_index.push_back(Array1d<int>(delta_t));
		dt_middle_index.push_back(Array1d<int>(delta_t));
		delta_ts.push_back(delta_t);
	}

	void save_subband_values(int idt, int _dt_middle_index, int _dt_middle_larger, int _dt_rest_index) {
		int idx = dt_middle_larger.size() - 1;
		dt_middle_larger[idx][idt] = _dt_middle_larger;
		dt_middle_index[idx][idt] = _dt_middle_index;
		dt_rest_index[idx][idt] = _dt_rest_index;
	}

	void copy_to_device() {
		for(int i = 0 ; i < dt_middle_larger.size(); i++) {
			dt_middle_larger[i].copy_to_device();
			dt_middle_index[i].copy_to_device();
			dt_rest_index[i].copy_to_device();
		}
	}
};



typedef struct _fdmt_t
{
	float fmin; // MHz
	float fmax; // MHz
	float df; // Channel offsets
	int order; // Order of the FMDT == log2(nf)
	int max_dt; // Maximum number of time integrations
	int nf; // Number of frequency bins
	int delta_t; // Initial delta_t for the initialisation
	int nbeams; // Number of beams
	int verbose; // 1 for verbose
	array4d_t states[2]; // iteration states
	int state_nbytes; // number of bytes in state
	int state_size; // number of elements in state
	FdmtIteration iterations[MAX_ITER]; // description of what happens for each iteration
} fdmt_t;

float dm_delay(const float f1, const float f2) ;

__host__ __device__ float squaref(const float f);

__host__ __device__ float isquaref(const float f);

__host__ __device__ float cff(float f1_start, float f1_end, float f2_start, float f2_end);

__host__ __device__ int calc_delta_t(const fdmt_t* fdmt, float f_start, float f_end);

int fdmt_create(fdmt_t* fdmt, float fmin, float fmax, int nf, int max_dt, int nbeams);

int fdmt_execute(fdmt_t* fdmt, fdmt_dtype* indata, fdmt_dtype* outdata);

/**
 * Does an FDMT iteration
 *
 * fdmt - instance
 * indata - 3D array, with dimensions [nf, n_d1, nt]
 * iteration_num - INteration number
 * output 3D array with dimensions [nf/2, n_d2, nt]
 */
int fdmt_iteration(const fdmt_t* fdmt,
		const int iteration_num,
		const array4d_t* indata,
		array4d_t* outdata);



int fdmt_initialise(const fdmt_t* fdmt, const array3d_t* indata, array4d_t* state);




#endif

