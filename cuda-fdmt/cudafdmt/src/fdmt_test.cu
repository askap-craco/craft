//
//  fdmt_test.c
//  fdmt
//
//  Created by Keith Bannister on 19/07/2016.
//  Copyright (c) 2016 Keith Bannister. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <iostream>
#include "fdmt.h"
#include "array.h"
#include "boxcar.h"
#include "CudaTimer.h"
#include "SigprocFile.h"
#include "rescale.h"


using namespace std;

void runtest_usage() {
	fprintf(stderr,
			"fdmt_test [options] infile outfile\n"
			"	-d Number of dispersion trials\n"
			"	-t Samples per block\n"
			"	-f Number of frequency channels\n"
			"	-b Number of beams\n"
			"	-x Maximum frequency (MHz)\n"
			"	-h Print this message\n"
	);
	exit(EXIT_FAILURE);
}

int main(int argc, char* argv[])
{
	printf("Test!");
	int nd = 512;
	int nt = 256;
	const int num_skip_blocks = 4;
	const int num_rescale_blocks = 2;
	float decay_timescale = 0.2; // Seconds?
	char ch;
	float thresh = 10.0;
	const char* out_filename = "fredda.cand";
	bool dump_data = false;
	while ((ch = getopt(argc, argv, "d:t:s:o:x:Dh")) != -1) {
		switch (ch) {
		case 'd':
			nd = atoi(optarg);
			break;
		case 't':
			nt = atoi(optarg);
			break;
		case 's':
			decay_timescale = atof(optarg);
			break;
		case 'o':
			out_filename = optarg;
			break;
		case 'x':
			thresh = atof(optarg);
			break;
		case 'D':
			dump_data = true;
			break;
		case '?':
		case 'h':
		default:
			runtest_usage();
		}
	}
	argc -= optind;
	argv += optind;

	if (argc == 0) {
		printf("Not enough arguments\n");
		exit(EXIT_FAILURE);
	}

	// Load sigproc file
	SigprocFile spf(argv[0]);
	CandidateSink sink(&spf, out_filename);
	cout << "spf tsamp " << spf.header_double("tsamp") << " nifs " << spf.header_int("nifs") << " fch1 " << spf.header_double("fch1")
							 << "foff " << spf.header_double("foff") << endl;
	int nbeams = spf.m_nifs;
	int nf = spf.m_nchans;
	size_t in_chunk_size = nbeams*nf*nt;

	// Create read buffer
	uint8_t* read_buf = (uint8_t*) malloc(sizeof(uint8_t) * in_chunk_size);
	array4d_t read_arr;
	read_arr.nw = 1;
	read_arr.nx = nt;
	read_arr.ny = nbeams;
	read_arr.nz = nf;

	assert(read_buf);

	array4d_t rescale_buf;
	rescale_buf.nw = nbeams;
	rescale_buf.nx = nf;
	rescale_buf.ny = 1;
	rescale_buf.nz = nt;
	array4d_malloc_hostonly(&rescale_buf);

	array4d_t out_buf;
	out_buf.nw = nbeams;
	out_buf.nx = 1;
	out_buf.ny = nd;
	out_buf.nz = nt;
	array4d_malloc_hostonly(&out_buf);

	// create rescaler
	rescale_t rescale;
	rescale.target_mean = 0.0;
	rescale.target_stdev = 1.0/sqrt((float) nf);
	rescale.decay_constant = 0.35 * decay_timescale / spf.m_tsamp; // This is how the_decimator.C does it, I think.
	rescale_allocate(&rescale, nbeams*nf);

	float foff =  (float) spf.m_foff;
	assert(foff < 0);
	float fmax = (float) spf.m_fch1 - foff; // The FDMT seems to want this to make sense of the world. Not sure why.
	float fmin = fmax + nf*foff;
	fdmt_t fdmt;
	printf("Creating FDMT fmin=%f fmax=%f nf=%d nd=%d nt=%d nbeams=%d\n", fmin, fmax, nf, nd, nt, nbeams);
	fdmt_create(&fdmt, fmin, fmax, nf, nd, nt, nbeams);
	spf.seek_sample(num_skip_blocks*nt);
	int blocknum = 0;

	while (spf.read_samples_uint8(nt, read_buf) == nt) {
		// File is in TBF order
		// Output needs to be BFT order
		// Do transpose and cast to float on the way through
		// TODO: Optimisation: cast to float and do rescaling in SIMD
		for(int t = 0; t < nt; ++t) {
			for (int b = 0; b < nbeams; ++b) {
				for (int f = 0; f < nf; ++f) {
					// NOTE: FDMT expects channel[0] at fmin
					// so invert the frequency axis if the frequency offset is negative
					int outf = f;
					if (spf.m_foff < 0) {
						outf = nf - f - 1;
					}
					int inidx = array4d_idx(&read_arr, 0, t, b, f);
					int outidx = array4d_idx(&rescale_buf, b, outf, 0, t);

					//printf("t=%d b=%d f=%d inidx=%d outidx=%d\n", t, b, f, inidx, outidx);
					// writes to inbuf
					size_t rs_idx = outf + nf*b;
					float v_rescale;
					v_rescale = rescale_update_decay_float_single(&rescale, rs_idx, (float) read_buf[inidx]);
					//float v_rescale = (((float) read_buf[inidx]) - 128.0f)*0.125;
					//v_rescale = (float) read_buf[inidx] - 128.;
					if (outf > 312) {
						v_rescale = 0;
					}
					rescale_buf.d[outidx] = v_rescale;
				}
			}
			rescale.sampnum += 1; // WARNING: Need to do this because we're calling rescale*single. THink harder about how to do this beter
		}

		char fbuf[1024];
		if (dump_data) {
			sprintf(fbuf, "inbuf_e%d.dat", blocknum);
			array4d_dump(&rescale_buf, fbuf);
			printf("Dumping input buffer %s\n", fbuf);
		}

		if (blocknum % num_rescale_blocks == 0) {
			rescale_update_scaleoffset(&rescale);
		}

		if (blocknum > num_rescale_blocks) {
			fdmt_execute(&fdmt, rescale_buf.d, out_buf.d);
			if (dump_data) {
				sprintf(fbuf, "fdmt_e%d.dat", blocknum);
				printf("Dumping fdmt buffer %s\n", fbuf);
				array4d_dump(&out_buf, fbuf);
			}
			boxcar_threshonly(&out_buf, thresh, sink);
		}

		blocknum++;
	}

	printf("FREDDA Finished\n");
	fdmt_print_timing(&fdmt);
}
int runtest(int argc, char* argv[])
{
	int nd = 512;
	int nt = 256;
	int nf = 336;
	int nbeams = 1;
	float fmax = 1440;
	char ch;
	while ((ch = getopt(argc, argv, "d:t:f:b:x:g:h")) != -1) {
		switch (ch) {
		case 'd':
			nd = atoi(optarg);
			break;
		case 't':
			nt = atoi(optarg);
			break;
		case 'f':
			nf = atoi(optarg);
			break;
		case 'b':
			nbeams = atoi(optarg);
			break;
		case 'x':
			fmax = atof(optarg);
			break;
		case '?':
		case 'h':
		default:
			runtest_usage();
		}
	}
	argc -= optind;
	argv += optind;

	float fmin = fmax - (float)nf;

	int blockin = nf*nt;
	int blockout = nd*nt;
	fdmt_dtype* din = (fdmt_dtype*) malloc(sizeof(fdmt_dtype)*blockin*nbeams);
	fdmt_dtype* din_tmp = (fdmt_dtype*) malloc(sizeof(fdmt_dtype)*blockin*nbeams);
	fdmt_dtype* dout = (fdmt_dtype*) malloc(sizeof(fdmt_dtype)*blockout*nbeams);
	printf("Starting! fmin=%f fmax=%f nbeams=%d nf=%d nd=%d nt=%d\n", fmin, fmax, nbeams, nf, nd, nt);

	if (argc != 2) {
		printf("Not enough arguments\n");
		exit(EXIT_FAILURE);
	}

	FILE* fin = fopen(argv[0], "r");
	if (fin == NULL) {
		perror("Could not open input file");
		exit(EXIT_FAILURE);
	}

	FILE* fout = fopen(argv[1], "w");
	if (fout == NULL) {
		perror("Could not open output file");
		exit(EXIT_FAILURE);
	}


	fdmt_t fdmt;
	fdmt_create(&fdmt, fmin, fmax, nf, nd, nt, nbeams);

	int nbox = 32;
	array4d_t boxout;
	boxout.nw = nbeams;
	boxout.nx = nd;
	boxout.ny = nt;
	boxout.nz = nbox;
	array4d_malloc(&boxout);

	// read input file until exhausted
	while (fread(din_tmp, sizeof(fdmt_dtype), blockin, fin) == blockin) {

		// File is in TF format. We need FT order.
		// Do the transpose
		for(int t = 0; t < nt; ++t) {
			for (int f = 0; f < nf; f++) {
				din[f*nt + t] = din_tmp[f + nf*t];
			}
		}
		// copy to all beams
		for(int b = 1; b < nbeams; b++) {
			int idx = b*blockin;
			//memcpy(&din[idx], din, blockin*sizeof(fdmt_dtype));
		}

		CudaTimer t;
		t.start();
		for(int i = 0; i < 1; i++) {
			fdmt_execute(&fdmt, din, dout);
		}

		boxcar_do(&fdmt.states[fdmt.curr_state_idx], &boxout);

		t.stop();
		cout << "FDMT Execute loop took " << t << endl;
		fwrite(dout, sizeof(fdmt_dtype), blockout, fout);
		cout << "Wrote " << blockout << " elements to outfile. First two are:" << dout[0] << dout[1] << endl;
	}
	fclose(fin);
	fclose(fout);
}
