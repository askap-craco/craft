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
#include <limits.h>
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
#include "CandidateList.h"

#include "rescale.h"


using namespace std;

void runtest_usage() {
	fprintf(stderr,
			"cudafdmt [options] [infile [infile[ ...]]\n"
			"   -d D - Number of dispersion trials. Negative D computes negative DMs\n"
			"   -t T - Samples per block\n"
			"   -s S - Decay timescale\n"
			"   -o FILE - Candidate filename\n"
			"   -x SN - threshold S/N\n"
			"   -D dump intermediate data to disk\n"
			"   -r R - Blocks per rescale update (0 for no rescaling)\n"
			"   -S S - Seek to this number of seconds before starting\n"
			"   -M M - Channel Mean flagging threshold (3 is OK)\n"
			"   -T T - Channel StdDev flagging threshold (3 is OK)\n"
			"   -K K - Channel Kurtosis threshold (0.8 is pretty good)\n"
			"   -G N - Channel flag channel growing (flags N channels either side of a bad channel)\n"
			"   -z Z - Zap times with 0 DM above threshold Z\n"
			"   -C C - Zap time/frequency cells with S/N above threshold C\n"
			"   -u   - Subtract DM0 time series from spectrum\n"
			"   -n ncand - Maximum mumber of candidates to write per block\n"
			"   -m mindm - Minimum DM to report candidates for (to ignore 0 DM junk)\n"
			"   -b maxbc - Maximum boxcar to create a candidate. Candidates with peaks above this boxcar are ignored\n"
			"   -g G - CUDA device\n"
			"   -N N - Maximum number of blocks to process before quitting\n"
			"   -h Print this message\n"
			"    Version: %s\n"
	, VERSION);
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
	if (copy) {
		array4d_copy_to_host(arr);
	}
	int nz = 0;
	int size = array4d_size(arr);
	for(int i = 0; i < size; i++) {
		if (arr->d[i] == 0.0) {
			nz += 1;
		}
	}

	printf("Dumping %s %s %d zeros\n", prefix, fbuf, nz);
	array4d_dump(arr, fbuf);
}

int main(int argc, char* argv[])
{
	int nd = 1024;
	int nt = 512;
	float seek_seconds = 0.0;
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
	float dm0_thresh = 1e9;
	float cell_thresh = 1e9;
	int flag_grow = 3;
	int max_ncand_per_block = 4096;
	int mindm = 0;
	int maxbc = 32;
	int max_nblocks = INT_MAX;
	bool subtract_dm0 = false;

	printf("Fredda version %s starting. Cmdline: ", VERSION);
	for (int c = 0; c < argc; ++c) {
		printf("%s ", argv[c]);
	}
	printf("\n");

	while ((ch = getopt(argc, argv, "d:t:s:o:x:r:S:Dg:M:T:K:G:C:n:m:b:z:N:uh")) != -1) {
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
			seek_seconds = atof(optarg);
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
		case 'G':
			flag_grow = atoi(optarg);
			break;
		case 'C':
			cell_thresh = atof(optarg);
			break;
		case 'n':
			max_ncand_per_block = atoi(optarg);
			break;
		case 'm':
			mindm = atoi(optarg);
			break;
		case 'u':
			subtract_dm0 = true;
			break;
		case 'b':
			maxbc = atoi(optarg);
			break;
		case 'z':
			dm0_thresh = atof(optarg);
			break;
		case 'N':
			max_nblocks = atoi(optarg);
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

	printf("\n");
	printf("Setting cuda device to %d\n", cuda_device);
	gpuErrchk( cudaSetDevice(cuda_device));

	CpuTimer tall;
	CudaTimer trescale;
	CudaTimer tboxcar;
	tall.start();

	// Load sigproc file
	SigprocFileSet source(argc, argv);
	bool negdm = (nd < 0);

	CandidateSink sink(&source, out_filename, negdm);
	cout << "spf tsamp " << source.tsamp()<< " nbeams " << source.nbeams() << " fch1 " << source.fch1() << " nchans "
			<< source.nchans() << "foff " << source.foff() << endl;
	int nbeams = source.nbeams();
	int nf = source.nchans();
	size_t in_chunk_size = nbeams*nf*nt;

	float foff =  (float) source.foff();
	assert(foff < 0);
	float fmax = (float) source.fch1() - foff; // The FDMT seems to want this offset to make sense of the world. Not sure why.
	float fmin = fmax + nf*foff;


	if (nd < 0) { // Flip the band to calculate negative DMs
		nd = -nd; // make nd positive -otherwise array sizes get confuddled
		// FDMT requres fmin < fmax
		// rescaling will invert the channels now that we've changed the sign of foff
		foff = -foff;
	}

	// Create read buffer
	uint8_t* read_buf = (uint8_t*) malloc(sizeof(uint8_t) * in_chunk_size);
	assert(read_buf);

	array4d_t read_arr;
	read_arr.nw = 1;
	read_arr.nx = nt;
	read_arr.ny = nbeams;
	read_arr.nz = nf;

	array4d_t rescale_buf;
	rescale_buf.nw = nbeams;
	rescale_buf.nx = nf;
	rescale_buf.ny = 1;
	rescale_buf.nz = nt;
	array4d_malloc(&rescale_buf, dump_data, true);

	array4d_t out_buf;
	out_buf.nw = nbeams;
	out_buf.nx = 1;
	out_buf.ny = nd;
	out_buf.nz = nt;
	array4d_malloc(&out_buf, dump_data, true);

	// create rescaler
	rescale_gpu_t rescale;
	rescale.interval_samps = nt;
	rescale.target_mean = 0.0;
	//rescale.target_stdev = 1.0/sqrt((float) nf);
	rescale.target_stdev = 1.0;
	rescale.decay_constant = 0.35 * decay_timescale / source.tsamp(); // This is how the_decimator.C does it, I think.
	rescale.mean_thresh = mean_thresh;
	rescale.std_thresh = std_thresh;
	rescale.kurt_thresh = kurt_thresh;
	rescale.flag_grow = flag_grow;
	rescale.dm0_thresh = dm0_thresh;
	rescale.cell_thresh = cell_thresh;
	// set guess of initial scale and offset to dm0 thresholding works
	printf("Rescaling to mean=%f stdev=%f decay constant=%f mean/std/kurtosis/dm0/Cell thresholds: %0.1f/%0.1f/%0.1f/%0.1f/%0.1f grow flags by %d channels\n",
			rescale.target_mean,rescale.target_stdev,
			rescale.decay_constant,
			rescale.mean_thresh, rescale.std_thresh, rescale.kurt_thresh,
			rescale.dm0_thresh, rescale.cell_thresh,
			rescale.flag_grow);
	//rescale_allocate(&rescale, nbeams*nf);
	rescale_allocate_gpu(&rescale, nbeams, nf, nt, true); // Need host memory allocated for rescale because we copy back to count flags
	if (num_rescale_blocks == 0) {
		rescale_set_scale_offset_gpu(&rescale, 1.0f, -128.0f); // Just pass it straight through without rescaling
	} else {
		rescale_set_scale_offset_gpu(&rescale, rescale.target_stdev/18.0, -128.0f); // uint8 stdev is 18 and mean +128.
	}

	fdmt_t fdmt;
	printf("Creating FDMT fmin=%f fmax=%f nf=%d nd=%d nt=%d nbeams=%d\n", fmin, fmax, nf, nd, nt, nbeams);
	fdmt_create(&fdmt, fmin, fmax, nf, nd, nt, nbeams, dump_data);
	assert(seek_seconds >= 0);
	int num_skip_blocks = seek_seconds / source.tsamp() / nt;
	printf("Seeking to start of data: block %d nsamples=%d time=%fs\n", num_skip_blocks, num_skip_blocks*nt, num_skip_blocks*nt*source.tsamp());
	printf("S/N Threshold %f Max ncand per block %d mindm %d \n", thresh, max_ncand_per_block, mindm);
	source.seek_sample(nt*num_skip_blocks);
	int blocknum = 0;
	int iblock = num_skip_blocks;
	unsigned long long total_candidates = 0;

	// make boxcar history
	array4d_t boxcar_history;
	boxcar_history.nw = 1;
	boxcar_history.nx = nbeams;
	boxcar_history.ny = nd;
	boxcar_history.nz = NBOX;
	array4d_malloc(&boxcar_history, dump_data, true);
	array4d_zero(&boxcar_history);

	// make boxcar output.
	// TODO: Only allocate on GPU if we'll be dumping it to dis.
	// Otherwise, we'll just use candidate lists and save on a bucketload of memory
	array4d_t boxcar_data;
	boxcar_data.nw = nbeams;
	boxcar_data.nx = nd;
	boxcar_data.ny = nt;
	boxcar_data.nz = NBOX;
	array4d_malloc(&boxcar_data, dump_data, dump_data);
	array4d_zero(&boxcar_data);

	CandidateList candidate_list(max_ncand_per_block);

	// measure bytes used
	size_t gpu_free_bytes, gpu_total_bytes;
	gpuErrchk(cudaMemGetInfo( &gpu_free_bytes, &gpu_total_bytes ));

	// add signal handler
	signal(SIGHUP, &handle_signal);
	signal(SIGINT, &handle_signal);
	int num_flagged_beam_chans = 0;
	int num_flagged_times = 0;
	int blocks_since_rescale_update = 0;

	while (source.read_samples_uint8(nt, read_buf) == nt) {
		if (stopped) {
			printf("Stopped due to signal received\n");
			break;
		}
		if (blocknum >= max_nblocks) {
			break;
		}

		// File is in TBF order
		// Output needs to be BFT order
		// Do transpose and cast to float on the way through using GPU
		// copy raw data to state. Here we're a little dodgey

		bool invert_freq = (foff < 0);
		uint8_t* read_buf_device = (uint8_t*) fdmt.states[0].d_device;
		fdmt.t_copy_in.start();
		gpuErrchk(cudaMemcpy(read_buf_device, read_buf, in_chunk_size*sizeof(uint8_t), cudaMemcpyHostToDevice));
		fdmt.t_copy_in.stop();
		trescale.start();
		rescale_update_and_transpose_float_gpu(rescale, rescale_buf, read_buf_device, invert_freq, subtract_dm0);
		trescale.stop();

		if (dump_data) {
			dumparr("inbuf", iblock, &rescale_buf);
		}

		// Count how many times were flagged
		assert(num_rescale_blocks >= 0);
		array4d_copy_to_host(&rescale.nsamps); // must do this before updaing scaleoffset, which resets nsamps to zero

		for(int i = 0; i < nf*nbeams; ++i) {
			int nsamps = (int)rescale.nsamps.d[i]; // nsamps is the number of unflagged samples from this block
			int nflagged = rescale.sampnum - nsamps;
			// rescale.sampnum is the total number of samples that has gone into the rescaler
			assert (nflagged >= 0);
			num_flagged_times += nflagged;
		}



		// do rescaling if required
		if (num_rescale_blocks > 0 && blocknum % num_rescale_blocks == 0) {
			rescale_update_scaleoffset_gpu(rescale);

			// Count how many  channels have been flagged for this whole block
			// by looking at how many channels have scale==0
			array4d_copy_to_host(&rescale.scale);
			for(int i = 0; i < nf*nbeams; ++i) {
				if (rescale.scale.d[i] == 0) {
					// that channel will stay flagged for num_rescale_blocks
					num_flagged_beam_chans += num_rescale_blocks;
				}
				// Count how many times have been flagged for this block
				// TODO: DANGER DANGER! This doesn't count flagged times if num_rescale_blocks = 0
				// This gave me a long headache at LAX when I set -s 1e30 stupidly.
				int nsamps = (int)rescale.nsamps.d[i];
				// nsamps is the number of unflagged samples in nt*num_rescale_blocks samples
				int nflagged = nt*num_rescale_blocks - nsamps;
				assert (nflagged >= 0);
				num_flagged_times += nflagged;
			}

			if (dump_data) {
				dumparr("mean", iblock, &rescale.mean);
				dumparr("std", iblock, &rescale.std);
				dumparr("kurt", iblock, &rescale.kurt);
				dumparr("nsamps", iblock, &rescale.nsamps);
				dumparr("dm0", iblock, &rescale.dm0);
				dumparr("dm0count", iblock, &rescale.dm0count);
				dumparr("dm0stats", iblock, &rescale.dm0stats);
			}
		}

		if (blocknum >= num_rescale_blocks) {
			/// Execute the FDMT
			fdmt_execute(&fdmt, rescale_buf.d_device, out_buf.d);
			if (dump_data) {
				dumparr("fdmt", iblock, &out_buf, false);
			}
			size_t sampno = iblock*nt;

			//total_candidates += boxcar_threshonly(&out_buf, sampno, thresh, max_ncand_per_block, mindm, sink);
			tboxcar.start();
			boxcar_do_gpu (
					&fdmt.ostate,
					&boxcar_data,
					&boxcar_history,
					thresh, max_ncand_per_block, mindm, maxbc, &candidate_list);
			tboxcar.stop();
			total_candidates += candidate_list.copy_to_sink(sink, sampno);

			if (dump_data) {
				dumparr("boxcar", iblock, &boxcar_data, true);
			}
		}

		blocknum++;
		iblock++;
	}

	float boxcar_ngops = nbeams*nt*nd*2*NBOX/1e9;

	float flagged_percent = ((float) num_flagged_beam_chans) / ((float) nf*nbeams*blocknum) * 100.0f;
	float dm0_flagged_percent = ((float) num_flagged_times) / ((float) blocknum*nbeams*nt*nf) * 100.0f;
	printf("FREDDA Finished\nFound %llu candidates \n", total_candidates);
	float data_nsecs = blocknum*nt*source.tsamp();
	tall.stop();
	cout << "Processed " << blocknum << " blocks = "<< blocknum*nt << " samples = " << data_nsecs << " seconds" << " at " << data_nsecs/tall.wall_total()<< "x real time"<< endl;
	cout << "Freq auto-flagged " << num_flagged_beam_chans << "/" << (nf*nbeams*blocknum) << " channels = " << flagged_percent << "%" << endl;
	cout << "DM0 auto-flagged " << num_flagged_times << "/" << (blocknum*nbeams*nt*nf) << " samples = " << dm0_flagged_percent << "%" << endl;
	cout << "FREDDA CPU "<< endl << tall << endl;
	cout << "Rescale "<< endl << trescale << endl;
	cout << "Boxcar "<< endl << tboxcar << endl;
	cout << "File reading " << endl << source.read_timer << endl;
	fdmt_print_timing(&fdmt);
	cout << "FDMT " << fdmt.nops/1e9
			<< " Gops/iteration ran at: " << fdmt.nops / (fdmt.t_iterations.get_average_time()/1e3)/1e9
			<< " GFLOPS" << endl;
	cout << "Boxcar " << boxcar_ngops
			<< " Gops/iteration. ran at: " << boxcar_ngops/(tboxcar.get_average_time()/1e3)
			<< " GFLOPS" << endl;
	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	cout << "Resources User: " << usage.ru_utime.tv_sec <<
			"s System:" << usage.ru_stime.tv_sec << "s MaxRSS:" << usage.ru_maxrss/1024/1024 << "MB" << endl;
	cout << "GPU Memory used " << (gpu_total_bytes - gpu_free_bytes)/1024/1024 << " of " << gpu_total_bytes /1024/124 << " MiB" << endl;
}

