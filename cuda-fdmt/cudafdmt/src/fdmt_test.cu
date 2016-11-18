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
#include <signal.h>
#if defined (ENABLE_OPENMP)
#include <omp.h>
#endif
#include <sys/time.h>
#include <sys/resource.h>
#include "fdmt.h"
#include "array.h"
#include "boxcar.h"
#include "CudaTimer.h"
#include "CpuTimer.h"
#include "DataSource.h"
#include "SigprocFile.h"
#include "SigprocFileSet.h"

#include "rescale.h"


using namespace std;

void runtest_usage() {
	fprintf(stderr,
			"cudafdmt [options] [infile [infile[ ...]]\n"
			"   -d Number of dispersion trials\n"
			"   -t Samples per block\n"
			"   -s Decay timescale\n"
			"   -o Candidate filename\n"
			"   -x threshold S/N\n"
			"   -D dump intermediate data to disk\n"
			"   -r Blocks per rescale update\n"
			"   -S Number of blocks to skip\n"
			"   -K Kurtosis threshold (0.8 is pretty good)\n"
			"   -M Mean offset threshold (3 is OK)\n"
			"   -T StdDev threshold (3 is OK)\n"
			"   -G Flag channel growing (flags N channels either side of a bad channel)\n"
			"   -g CUDA device\n"
			"   -h Print this message\n"
	);
	exit(EXIT_FAILURE);
}

volatile bool stopped;

//typedef void (*sig_t) (int);

void handle_signal(int signal)
{
	stopped = true;
}

void dumparr(const char* prefix, const int blocknum, array4d_t* arr, bool copy=true)
{
	char fbuf[1024];
	sprintf(fbuf, "%s_e%d.dat", prefix, blocknum);
	printf("Dumping %s %s\n", prefix, fbuf);
	if (copy) {
		array4d_copy_to_host(arr);
	}
	array4d_dump(arr, fbuf);
}

int main(int argc, char* argv[])
{
	int nd = 512;
	int nt = 256;
	int num_skip_blocks = 4;
	int num_rescale_blocks = 2;
	float decay_timescale = 0.2; // Seconds?
	char ch;
	float thresh = 10.0;
	const char* out_filename = "fredda.cand";
	bool dump_data = false;
	int cuda_device = 0;
	float kurt_thresh = 1e9;
	float std_thresh = 1e9;
	float mean_thresh = 1e9;
	int flag_grow = 3;
	while ((ch = getopt(argc, argv, "d:t:s:o:x:r:S:Dg:M:T:K:G:h")) != -1) {
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
		case 'r':
			num_rescale_blocks = atoi(optarg);
			break;
		case 'S':
			num_skip_blocks = atoi(optarg);
			break;
		case 'g':
			cuda_device = atoi(optarg);
			break;
		case 'K':
			kurt_thresh = atof(optarg);
			break;
		case 'T':
			std_thresh = atof(optarg);
			break;
		case 'M':
			mean_thresh = atof(optarg);
			break;
		case 'F':
			flag_grow = atoi(optarg);
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
		printf("Not enough arguments: %d\n");
		exit(EXIT_FAILURE);
	}

	printf("Setting cuda device to %d\n", cuda_device);
	gpuErrchk( cudaSetDevice(cuda_device));

	CpuTimer tall;
	CpuTimer trescale;
	CpuTimer tboxcar;
	tall.start();

	// Load sigproc file
	SigprocFileSet source(argc, argv);

	CandidateSink sink(&source, out_filename);
	cout << "spf tsamp " << source.tsamp()<< " nbeams " << source.nbeams() << " fch1 " << source.fch1() << " nchans "
			<< source.nchans() << "foff " << source.foff() << endl;
	int nbeams = source.nbeams();
	int nf = source.nchans();
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
	array4d_malloc(&rescale_buf); // Can do GPU only maybe??

	array4d_t out_buf;
	out_buf.nw = nbeams;
	out_buf.nx = 1;
	out_buf.ny = nd;
	out_buf.nz = nt;
	array4d_malloc_hostonly(&out_buf);

	// create rescaler
	rescale_gpu_t rescale;
	rescale.interval_samps = nt;
	rescale.target_mean = 0.0;
	rescale.target_stdev = 1.0/sqrt((float) nf);
	rescale.decay_constant = 0.35 * decay_timescale / source.tsamp(); // This is how the_decimator.C does it, I think.
	rescale.mean_thresh = mean_thresh;
	rescale.std_thresh = std_thresh;
	rescale.kurt_thresh = kurt_thresh;
	rescale.flag_grow = flag_grow;
	printf("Rescaling to mean=%f stdev=%f decay constant=%f mean/std/kurtosis thresholds: %f/%f/%f grow flags by %d channels\n",rescale.target_mean,rescale.target_stdev,
			rescale.decay_constant, rescale.mean_thresh, rescale.std_thresh, rescale.kurt_thresh,rescale.flag_grow);
	//rescale_allocate(&rescale, nbeams*nf);
	rescale_allocate_gpu(&rescale, nbeams, nf);
	rescale_set_scale_offset_gpu(&rescale, 1.0, 0.0); // Set initial scale and offset

	float foff =  (float) source.foff();
	assert(foff < 0);
	float fmax = (float) source.fch1() - foff; // The FDMT seems to want this to make sense of the world. Not sure why.
	float fmin = fmax + nf*foff;
	fdmt_t fdmt;
	printf("Creating FDMT fmin=%f fmax=%f nf=%d nd=%d nt=%d nbeams=%d\n", fmin, fmax, nf, nd, nt, nbeams);
	fdmt_create(&fdmt, fmin, fmax, nf, nd, nt, nbeams);
	printf("Seeking to start of data: nblocks=%d nsamples=%d time=%fs\n", num_skip_blocks, num_skip_blocks*nt, num_skip_blocks*nt*source.tsamp());
	source.seek_sample(num_skip_blocks*nt);
	int blocknum = 0;

	// add signal handler
	signal(SIGHUP, &handle_signal);
	signal(SIGINT, &handle_signal);

	while (source.read_samples_uint8(nt, read_buf) == nt) {
		if (stopped) {
			printf("Stopped due to signal received\n");
			break;
		}
		//size_t nt2 = fin.read_samples_uint8(nt, read_buf2);
		//assert(nt2 = nt);
		// File is in TBF order
		// Output needs to be BFT order
		// Do transpose and cast to float on the way through
		// TODO: Optimisation: cast to float and do rescaling in SIMD

		// copy raw data to state. Here we're a little dodgey

		bool invert_freq = (foff < 0);
		array4d_t* inarr = &fdmt.states[1];
		uint8_t* read_buf_device = (uint8_t*) fdmt.states[0].d_device;
		fdmt.t_copy_in.start();
		gpuErrchk(cudaMemcpy(read_buf_device, read_buf, in_chunk_size*sizeof(uint8_t), cudaMemcpyHostToDevice));
		fdmt.t_copy_in.stop();
		trescale.start();
		//rescale_update_and_transpose_float(rescale, read_arr, rescale_buf,read_buf, invert_freq);
		rescale_update_and_transpose_float_gpu(rescale, rescale_buf, read_buf_device, invert_freq);
		trescale.stop();

		if (dump_data) {
			dumparr("inbuf", blocknum, &rescale_buf);
		}

		assert(num_rescale_blocks >= 0);

		if (num_rescale_blocks > 0 && blocknum % num_rescale_blocks == 0) {
			rescale_update_scaleoffset_gpu(rescale);
			if (dump_data) {
				dumparr("mean", blocknum, &rescale.mean);
				dumparr("std", blocknum, &rescale.std);
				dumparr("kurt", blocknum, &rescale.kurt);
			}
		}

		if (blocknum >= num_rescale_blocks) {
			fdmt_execute(&fdmt, rescale_buf.d_device, out_buf.d);
			if (dump_data) {
				dumparr("fdmt", blocknum, &out_buf, false);
			}
			tboxcar.start();
			boxcar_threshonly(&out_buf, thresh, sink);
			tboxcar.stop();

		}

		blocknum++;
	}

	printf("FREDDA Finished\n");
	tall.stop();
	cout << "Processed " << blocknum*nt << " samples = " << blocknum*nt*source.tsamp() << " seconds" << endl;
	cout << "FREDDA CPU "<< tall << endl;
	cout << "Rescale CPU "<< trescale << endl;
	cout << "Boxcar CPU "<< tboxcar << endl;
	cout << "File reading " << source.read_timer << endl;
	fdmt_print_timing(&fdmt);

	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	cout << "Resources User: " << usage.ru_utime.tv_sec <<
			"s System:" << usage.ru_stime.tv_sec << "s MaxRSS:" << usage.ru_maxrss/1024/1024 << "MB" << endl;
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
