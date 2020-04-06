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
#include "CudaTimer.h"
//#define DUMP_STATE 1

#define MAX_ITER 16

using namespace std;

typedef struct _FdmtSum {
	int id1;
	int id2;
	int toff;
	int inchan1;
	int inchan2;
	int iterno;
} FdmtDmConfig;

class FdmtSubbandConfig {
public:
	int delta_t;
	vector<FdmtDmConfig> dmconfig;
};

class FdmtIteration
{
public:

	vector<Array2d<int>* > dt_data;
	vector<int> delta_ts;
	coord4_t state_shape;
	vector<FdmtSubbandConfig*> subband_config;
	int m_iterno;
	int m_nbeam;
	int m_nf;
	int m_ndt;
	int m_nt;

	FdmtIteration(int iterno, int nbeam, int nf, int ndt, int nt) :
		m_iterno(iterno),
		m_nbeam(nbeam), m_nf(nf), m_ndt(ndt), m_nt(nt) {
		state_shape.w = nbeam;
		state_shape.x = nf;
		state_shape.y = ndt;
		state_shape.z = nt;
	}


	__host__ void add_subband(int delta_t) {
		dt_data.push_back(new Array2d<int>(delta_t, 4));
		delta_ts.push_back(delta_t);
		FdmtSubbandConfig* config = new FdmtSubbandConfig();
		config->delta_t = delta_t;
		subband_config.push_back(config);
	}

	__host__ void save_subband_values(int idt, int src1_offset, int src2_offset, int out_offset, int mint,
			int id1, int id2, int toff) {
		Array2d<int>* dts = dt_data.back();
		dts->set_host(idt, 0, src1_offset);
		dts->set_host(idt, 1, src2_offset);
		dts->set_host(idt, 2, out_offset);
		dts->set_host(idt, 3, mint);

		FdmtSubbandConfig* conf = subband_config.back();
		FdmtDmConfig dmconfig;
		dmconfig.id1 = id1;
		dmconfig.id2 = id2;
		dmconfig.toff = toff;
		conf->dmconfig.push_back(dmconfig);
	}

	__host__ void copy_to_device() {
		for(int i = 0 ; i < dt_data.size(); i++) {
			Array2d<int>* dts = dt_data.at(i);
			dts->copy_to_device();
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
	int nt; // Number of integrations per block
	int nbox; // Number of boxcar trials to calculate weights for
	int delta_t; // Initial delta_t for the initialisation
	int nbeams; // Number of beams
	int nbeams_alloc; // Number of beams to allocate memory for, and process simultaneously
	int verbose; // 1 for verbose
	int curr_state_idx; // index of the current state in the states[] array
	array4d_t states[2]; // iteration states
	array4d_t ostate; // Special state for output delay and sum
	array4d_t weights; // weights to apply to dispersion axis
	int state_nbytes; // number of bytes in state
	int state_size; // number of elements in state
	int execute_count; //  number of times execute() has been called
	float _df_top; // Width of the top subband of an FDMT iteration - used during fdmt_create() to handle non-power-of-two FDMTs
	float _df_bot; // Width of all other subbands of an FDMT iteration - used during fdmt_create() to handle non-power-of-two FDMTs
	int _ndt_top; // Number of valid dispersion measures in the top band
	//FdmtIteration* iterations[MAX_ITER]; // description of what happens for each iteration
        long long int nops; // Number of floating point operations we do for all iterations for a single beam
	vector<FdmtIteration* > iterations;
	bool dump_data;
	CudaTimer t_init;
	CudaTimer t_iterations;
	CudaTimer t_copy_in;
	CudaTimer t_update_ostate;
	CudaTimer t_copy_back;

} fdmt_t;

float dm_delay(const float f1, const float f2) ;

__host__ __device__ float squaref(const float f);

__host__ __device__ float isquaref(const float f);

__host__ __device__ float cff(float f1_start, float f1_end, float f2_start, float f2_end);

__host__ __device__ int calc_delta_t(const fdmt_t* fdmt, float f_start, float f_end);

int fdmt_create(fdmt_t* fdmt, float fmin, float fmax, int nf, int max_dt, int nt,  int nbeams, int nbeams_alloc, bool dump_data);

int fdmt_execute(fdmt_t* fdmt, fdmt_dtype* indata, fdmt_dtype* outdata);

int fdmt_calculate_weights(fdmt_t* fdmt);


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

void fdmt_print_timing(fdmt_t* fdmt);

#endif

