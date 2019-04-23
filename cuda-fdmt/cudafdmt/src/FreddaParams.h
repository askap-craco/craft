/*
 * FreddaParams.h
 *
 *  Created on: 15 Apr 2019
 *      Author: ban115
 */

#ifndef FREDDAPARAMS_H_
#define FREDDAPARAMS_H_
#include "DataSource.h"

class FreddaParams {
public:

	int nd = 1024;
	int nt = 512;
	float seek_seconds = 0.0;
	int num_rescale_blocks = 2;
	float decay_timescale = 1.0; // Seconds?
	char ch;
	float thresh = 10.0;
	const char* out_filename = "fredda.cand";
	const char* flag_file = NULL;
	bool dump_data = false;
	bool do_dump_rescaler = false;
	int cuda_device = 0;
	float kurt_thresh = INFINITY;
	float std_thresh = INFINITY;
	float mean_thresh = INFINITY;
	float dm0_thresh = INFINITY;
	float cell_thresh = INFINITY;
	int flag_grow = 3;
	int max_ncand_per_block = 4096;
	int mindm = 0;
	int maxbc = 32;
	int max_nblocks = INT_MAX;
	int nbeams_alloc = -1;
	bool subtract_dm0 = false;
	bool polsum = false;
	char udp_host[128];
	short udp_port = -1;
	int export_dada_key = -1;

	DataSource* source;
	double out_freq;
	double out_foff;
	int nbeams_per_antenna; // number of beams including polarisations
	int nbeams_in_total;
	int npols_in ;
	int nbeams_out;
	int npols_out;
	float nbeams_summed ;
	int nf;
	int nbits ;
	float foff;
	float fmax; // The FDMT seems to want this offset to make sense of the world. Not sure why.
	float fmin;
	int argc;
	char** argv;
	FreddaParams();
	virtual ~FreddaParams();
	void set_source(DataSource& source);
	void parse(int argc, char* argv[]);
	void to_dada(char header_buf[]);
};


#endif /* FREDDAPARAMS_H_ */
