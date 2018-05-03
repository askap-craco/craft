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
#include <float.h>
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
#include "DataSource.h"
#include "DadaSource.h"
#include "CandidateList.h"
#include "InvalidSourceFormat.h"
#include "Rescaler.h"


#include "rescale.h"


using namespace std;

void runtest_usage() {
	fprintf(stderr,
			"cudafdmt [options] [infile [infile[ ...]]\n"
			"   -d D - Number of dispersion trials. Negative D computes negative DMs\n"
			"   -t T - Samples per block\n"
			"   -s S - Decay timescale\n"
			"   -o FILE - Candidate filename\n"
			"   -U host:port - UDP host:port to send candidates to\n"
			"   -x SN - threshold S/N\n"
			"   -D dump intermediate data to disk\n"
			"   -B b - Process b beams simultaneously to save memory\n"
			"   -r R - Blocks per rescale update (0 for no rescaling)\n"
			"   -S S - Seek to this number of seconds before starting\n"
			"   -M M - Channel Mean flagging threshold (3 is OK)\n"
			"   -T T - Channel StdDev flagging threshold (3 is OK)\n"
			"   -K K - Channel Kurtosis threshold (0.8 is pretty good)\n"
			"   -G N - Channel flag channel growing (flags N channels either side of a bad channel)\n"
			"   -z Z - Zap times with 0 DM above threshold Z\n"
			"   -C C - Zap time/frequency cells with S/N above threshold C\n"
			"   -u   - Subtract DM0 time series from spectrum\n"
			"   -p   - Sum polarisations\n"
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
	bzero(udp_host, 128);
	short udp_port = -1;

	printf("Fredda version %s starting. Cmdline: ", VERSION);
	for (int c = 0; c < argc; ++c) {
		printf("%s ", argv[c]);
	}
	printf("\n");

	while ((ch = getopt(argc, argv, "d:t:s:o:x:r:S:B:Dg:M:T:U:K:G:C:n:m:b:z:N:uhp")) != -1) {
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
		case 'p':
			polsum = true;
			break;
		case 'B':
			nbeams_alloc = atoi(optarg);
			break;
		case 'U':
		{
			char* colon = strchr(optarg, ':');
			if (colon == NULL) {
				printf("Invalid hostport\n");
				exit(EXIT_FAILURE);
			}
			memcpy(udp_host, optarg, colon-optarg);
			udp_port = atoi(colon+1);
		}
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
	CpuTimer tproc;
	CudaTimer trescale;
	CudaTimer tboxcar;
	tall.start();

	DataSource* source = NULL;
	DadaSource* dada_source = NULL; // for debugging
	try {
		// load sigproc file
		SigprocFileSet* fs_source = new SigprocFileSet(nt, argc, argv);
		source = fs_source;
	} catch (InvalidSourceFormat& e) {
		try {
			int key;
			sscanf(argv[0], "%x",&key);
			dada_source = new DadaSource(nt, key, true);
			source = dada_source;
		} catch (InvalidSourceFormat& e) {
			printf("No valid inputs\n");
			exit(EXIT_FAILURE);
		}
	}
	assert(source != NULL);

	bool negdm = (nd < 0);
	CandidateSink sink(source, out_filename, negdm, udp_host, udp_port);
	cout << "spf tsamp " << source->tsamp()<< " nbeams " << source->nbeams()
			<< " npols "<< source->npols() << " fch1 " << source->fch1() << " nchans "
			<< source->nchans() << " foff " << source->foff() << endl;
	int nbeams_in = source->nbeams()*source->npols();
	int npols_in = source->npols();
	int nbeams_out;
	if (polsum) {
		nbeams_out = source->nbeams();
		assert(nbeams_in %2 == 0);
	} else {
		nbeams_out = nbeams_in;
	}
	float nbeams_summed = (float(nbeams_in)/float(nbeams_out));
	int nf = source->nchans();
	int nbits = source->nbits();
	size_t in_buffer_bytes = nbeams_in*nf*nt*nbits/8;
	void* in_buffer_device;
	gpuErrchk( cudaMalloc((void**) &in_buffer_device, in_buffer_bytes ));

	float foff =  (float) source->foff();
	assert(foff < 0);
	float fmax = (float) source->fch1() - foff; // The FDMT seems to want this offset to make sense of the world. Not sure why.
	float fmin = fmax + nf*foff;


	if (nd < 0) { // Flip the band to calculate negative DMs
		nd = -nd; // make nd positive -otherwise array sizes get confuddled
		// FDMT requres fmin < fmax
		// rescaling will invert the channels now that we've changed the sign of foff
		foff = -foff;
	}

	array4d_t rescale_buf;
	rescale_buf.nw = nbeams_out;
	rescale_buf.nx = nf;
	rescale_buf.ny = 1;
	rescale_buf.nz = nt;
	array4d_malloc(&rescale_buf, dump_data, true);

	array4d_t out_buf;
	out_buf.nw = nbeams_out;
	out_buf.nx = 1;
	out_buf.ny = nd;
	out_buf.nz = nt;
	array4d_malloc(&out_buf, dump_data, true);

	// create rescaler
	RescaleOptions rescale = {};
	rescale.interval_samps = nt;
	rescale.target_mean = 0.0;
	rescale.target_stdev = 1.0/sqrt(nbeams_summed);
	rescale.decay_constant = 0.35 * decay_timescale / source->tsamp(); // This is how the_decimator.C does it, I think.
	rescale.mean_thresh = mean_thresh;
	rescale.std_thresh = std_thresh;
	rescale.kurt_thresh = kurt_thresh;
	rescale.flag_grow = flag_grow;
	rescale.dm0_thresh = dm0_thresh;
	rescale.cell_thresh = cell_thresh;
	rescale.invert_freq = (foff < 0);
	rescale.subtract_dm0 = subtract_dm0;
	rescale.nt = nt;
	rescale.nf = nf;
	rescale.nbeams = nbeams_in;
	rescale.npols = npols_in;
	rescale.polsum = polsum;
	rescale.nbits = source->nbits();
	// set guess of initial scale and offset to dm0 thresholding works
	printf("Rescaling to mean=%f stdev=%f decay constant=%f mean/std/kurtosis/dm0/Cell thresholds: %0.1f/%0.1f/%0.1f/%0.1f/%0.1f grow flags by %d channels\n",
			rescale.target_mean,rescale.target_stdev,
			rescale.decay_constant,
			rescale.mean_thresh, rescale.std_thresh, rescale.kurt_thresh,
			rescale.dm0_thresh, rescale.cell_thresh,
			rescale.flag_grow);
	Rescaler* rescaler = new Rescaler(rescale);
	if (num_rescale_blocks == 0) {
		rescaler->set_scaleoffset(1.0f, 0.0f); // Just pass it straight through without rescaling
	}

	// Create fdmt
	fdmt_t fdmt;
	printf("Creating FDMT fmin=%f fmax=%f nf=%d nd=%d nt=%d nbeams=%d nbeams_alloc=%d\n",
			fmin, fmax, nf, nd, nt, nbeams_out, nbeams_alloc);
	fdmt_create(&fdmt, fmin, fmax, nf, nd, nt, nbeams_out, nbeams_alloc, dump_data);
	assert(seek_seconds >= 0);
	int num_skip_blocks = seek_seconds / source->tsamp() / nt;
	printf("Seeking to start of data: block %d nsamples=%d time=%fs\n", num_skip_blocks, num_skip_blocks*nt, num_skip_blocks*nt*source->tsamp());
	printf("S/N Threshold %f Max ncand per block %d mindm %d \n", thresh, max_ncand_per_block, mindm);
	source->seek_sample(nt*num_skip_blocks);
	int blocknum = 0;
	int iblock = num_skip_blocks;
	unsigned long long total_candidates = 0;
	unsigned long long num_candidate_overflow_blocks = 0;
	// make boxcar history
	array4d_t boxcar_history;
	boxcar_history.nw = 1;
	boxcar_history.nx = nbeams_out;
	boxcar_history.ny = nd;
	boxcar_history.nz = NBOX;
	array4d_malloc(&boxcar_history, dump_data, true);
	array4d_zero(&boxcar_history);
	// make boxcar discards
	array4d_t boxcar_discards;
	boxcar_discards.nw = 1;
	boxcar_discards.nx = 1;
	boxcar_discards.ny = nbeams_out;
	boxcar_discards.nz = nd;
	array4d_malloc(&boxcar_discards, true, true);
	array4d_cuda_memset(&boxcar_discards, 0);

	// make boxcar output.
	// TODO: Only allocate on GPU if we'll be dumping it to dis.
	// Otherwise, we'll just use candidate lists and save on a bucketload of memory
	array4d_t boxcar_data;
	boxcar_data.nw = nbeams_out;
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

	void* read_buf;

	while (source->read_samples(&read_buf) == nt) {
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
		// copy raw data to FDMT state. Here we're a little dodgey, but why allocate memory anyway?
		//uint8_t* read_buf_device = (uint8_t*) fdmt.states[0].d_device;
		fdmt.t_copy_in.start();
		gpuErrchk(cudaMemcpy(in_buffer_device, read_buf, in_buffer_bytes*sizeof(uint8_t), cudaMemcpyHostToDevice));
		fdmt.t_copy_in.stop();
		tproc.start();
		trescale.start();
		if (blocknum == 0) { // if first block rescale and update with no
			// flagging so we can work out roughly what the scales are
			rescaler->update_and_transpose(rescale_buf, in_buffer_device, rescaler->noflag_options);

			// update scale and offset
			rescaler->update_scaleoffset(rescaler->noflag_options);
		}

		// this time we rescale with the flagging turned on
		rescaler->update_and_transpose(rescale_buf, in_buffer_device, rescaler->options);
		trescale.stop();

		if (dump_data) {
			dumparr("inbuf", iblock, &rescale_buf);
		}

		// Count how many times were flagged
		assert(num_rescale_blocks >= 0);
		array4d_copy_to_host(&rescaler->nsamps); // must do this before updaing scaleoffset, which resets nsamps to zero

		for(int i = 0; i < nf*nbeams_in; ++i) {
			int nsamps = (int)rescaler->nsamps.d[i]; // nsamps is the number of unflagged samples from this block
			int nflagged = rescaler->sampnum - nsamps;
			// rescale.sampnum is the total number of samples that has gone into the rescaler
			assert (nflagged >= 0);
			num_flagged_times += nflagged;
		}

		// do rescaling if required
		if (num_rescale_blocks > 0 && blocknum % num_rescale_blocks == 0) {
			rescaler->update_scaleoffset(rescaler->options);

			// Count how many  channels have been flagged for this whole block
			// by looking at how many channels have scale==0
			array4d_copy_to_host(&rescaler->scale);
			for(int i = 0; i < nf*nbeams_in; ++i) {
				if (rescaler->scale.d[i] == 0) {
					// that channel will stay flagged for num_rescale_blocks
					num_flagged_beam_chans += num_rescale_blocks;
				}
				// Count how many times have been flagged for this block
				// TODO: DANGER DANGER! This doesn't count flagged times if num_rescale_blocks = 0
				// This gave me a long headache at LAX when I set -s 1e30 stupidly.
				int nsamps = (int)rescaler->nsamps.d[i];
				// nsamps is the number of unflagged samples in nt*num_rescale_blocks samples
				int nflagged = nt*num_rescale_blocks - nsamps;
				assert (nflagged >= 0);
				num_flagged_times += nflagged;
			}

			if (dump_data) {
				dumparr("mean", iblock, &rescaler->mean);
				dumparr("std", iblock, &rescaler->std);
				dumparr("kurt", iblock, &rescaler->kurt);
				dumparr("nsamps", iblock, &rescaler->nsamps);
				dumparr("dm0", iblock, &rescaler->dm0);
				dumparr("dm0count", iblock, &rescaler->dm0count);
				dumparr("dm0stats", iblock, &rescaler->dm0stats);
				dumparr("scale", iblock, &rescaler->scale);
				dumparr("offset", iblock, &rescaler->offset);
			}
		}

		if (blocknum >= num_rescale_blocks) {
			/// Execute the FDMT
			fdmt_execute(&fdmt, rescale_buf.d_device, out_buf.d);
			if (dump_data) {
				dumparr("fdmt", iblock, &out_buf, false);
				dumparr("ostate", iblock, & fdmt.ostate, true);
			}
			size_t sampno = iblock*nt;

			//total_candidates += boxcar_threshonly(&out_buf, sampno, thresh, max_ncand_per_block, mindm, sink);
			tboxcar.start();
			boxcar_do_gpu (
					&fdmt.ostate,
					&boxcar_data,
					&boxcar_history,
					&boxcar_discards,
					thresh, max_ncand_per_block, mindm, maxbc, &candidate_list);
			tboxcar.stop();
			int ncand = candidate_list.copy_to_sink(sink, sampno);
			if (ncand >= max_ncand_per_block - 1) {
				num_candidate_overflow_blocks++;
			}
			total_candidates += ncand;
			if (dump_data) {
				dumparr("boxcar", iblock, &boxcar_data, true);
			}
		}
		tproc.stop();

		blocknum++;
		iblock++;
	}

	tall.stop();

	// calculate array discards
	array4d_copy_to_host(&boxcar_discards);
	int total_discards = 0;
	for (int i = 0; i < array4d_size(&boxcar_discards); ++i) {
		total_discards += (int)boxcar_discards.d[i];
	}


	double boxcar_ngops = (double)nbeams_out*(double)nt*(double)nd*2.0*(double)NBOX/1e9;
	double data_nsecs = blocknum*nt*source->tsamp();

	double flagged_percent = ((double) num_flagged_beam_chans) / ((double) nf*nbeams_in*blocknum) * 100.0;
	double dm0_flagged_percent = ((double) num_flagged_times) / ((double) blocknum*nbeams_in*nt*nf) * 100.0;
	cout << " FREDDA Finished" << endl;
	cout << "Found " << total_candidates << " candidates" << endl;
	cout << "Discarded " << total_discards << " candidates for being too wide."<< endl;
	cout << num_candidate_overflow_blocks << " blocks overflowed the candidate buffer"<<endl;
	cout << "Processed " << blocknum << " blocks = "<< blocknum*nt << " samples = " << data_nsecs << " seconds" << " at " << data_nsecs/tall.wall_total()<< "x real time"<< endl;
	cout << "Freq auto-flagged " << num_flagged_beam_chans << "/" << (nf*nbeams_in*blocknum) << " channels = " << flagged_percent << "%" << endl;
	cout << "DM0 auto-flagged " << num_flagged_times << "/" << (blocknum*nbeams_in*nt*nf) << " samples = " << dm0_flagged_percent << "%" << endl;
	cout << "File reading " << endl << source->m_read_timer << endl;
	cout << "FREDDA Total "<< endl << tall << endl;
	cout << "FREDDA Procesing "<< endl << tproc << endl;
	cout << "Rescale "<< endl << trescale << endl;
	fdmt_print_timing(&fdmt);
	cout << "Boxcar "<< endl << tboxcar << endl;
	cout << "FDMT " << ((double)fdmt.nops)/1e9
			<< " Gops/iteration ran at: " << ((double)fdmt.nops) / (fdmt.t_iterations.get_average_time()/1e3)/1e9
			<< " GFLOPS" << endl;
	cout << "Boxcar " << boxcar_ngops
			<< " Gops/iteration. ran at: " << boxcar_ngops/(tboxcar.get_average_time()/1e3)
			<< " GFLOPS" << endl;
	struct rusage usage;
	getrusage(RUSAGE_SELF, &usage);
	cout << "Resources User: " << usage.ru_utime.tv_sec <<
			"s System:" << usage.ru_stime.tv_sec << "s MaxRSS:" << usage.ru_maxrss/1024/1024 << "MB" << endl;
	cout << "GPU Memory used " << (gpu_total_bytes - gpu_free_bytes)/1024/1024 << " of " << gpu_total_bytes /1024/124 << " MiB" << endl;
	delete source;
}

