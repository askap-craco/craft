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
#include <omp.h>
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
#include "DadaSet.h"
#include "FilDirSet.h"
#include "CandidateList.h"
#include "InvalidSourceFormat.h"
#include "Rescaler.h"
#include "rescale.h"
#include "DadaSink.h"
#include "FreddaParams.h"


using namespace std;


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

	//printf("Dumping %s %s %d zeros\n", prefix, fbuf, nz);
	array4d_dump(arr, fbuf);
}

void dump_rescaler(int iblock, Rescaler* rescaler)
{
	dumparr("mean", iblock, &rescaler->mean);
	dumparr("std", iblock, &rescaler->std);
	dumparr("kurt", iblock, &rescaler->kurt);
	dumparr("nsamps", iblock, &rescaler->nsamps);
	dumparr("dm0", iblock, &rescaler->dm0);
	dumparr("dm0count", iblock, &rescaler->dm0count);
	dumparr("dm0stats", iblock, &rescaler->dm0stats);
	dumparr("scale", iblock, &rescaler->scale);
	dumparr("offset", iblock, &rescaler->offset);
	dumparr("decay_offset", iblock, &rescaler->decay_offset);
}

int main(int argc, char* argv[])
{
	FreddaParams params; // this is new - I haven't finished refactoring everything to take advantage of FreddaParams
	params.parse(argc, argv);
	int nt = params.nt;
	int nd = params.nd;
	int nf = params.nf;
	printf("\n");
	printf("Setting cuda device to %d\n", params.cuda_device);
	gpuErrchk( cudaSetDevice(params.cuda_device));

	CpuTimer tall;
	CpuTimer tproc;
	CudaTimer trescale;
	CudaTimer tboxcar;
	tall.start();

	DataSource* source = NULL;
	DadaSet* dada_source = NULL; // for debugging
	try {
		// load sigproc file
		SigprocFileSet* fs_source = new SigprocFileSet(nt, argc, argv);
		source = fs_source;
	} catch (InvalidSourceFormat& e) {
		try {
			dada_source = new DadaSet(nt, argc, argv);
			source = dada_source;
		} catch (InvalidSourceFormat& e) {
			try {
				source = new FilDirSet(nt, argc, argv);
			} catch (InvalidSourceFormat& e) {
				printf("No valid inputs\n");
				exit(EXIT_FAILURE);
			}
		}
	}
	assert(source != NULL);
	int num_skip_blocks = params.seek_seconds / source->tsamp() / nt;
	printf("Seeking to start of data: block %d nsamples=%d time=%fs\n", num_skip_blocks, num_skip_blocks*nt, num_skip_blocks*nt*source->tsamp());
	if (num_skip_blocks > 0) {
		source->seek_sample(nt*num_skip_blocks);
	}
	params.set_source(*source);
	bool negdm = (params.nd < 0);
	CandidateSink sink(source, params.out_filename, negdm, params.udp_host, params.udp_port);
	cout << "spf tsamp " << source->tsamp()<< " ants " << source->nants() << " nbeams " << source->nbeams()
			<< " npols "<< source->npols() << " fch1 " << source->fch1() << " nchans "
			<< source->nchans() << " foff " << source->foff() << endl;

	printf("S/N Threshold %f Max ncand per block %d mindm %d \n", params.thresh, params.max_ncand_per_block, params.mindm);
	//rescale input buffer
	size_t in_buffer_bytes_per_ant = params.nbeams_per_antenna*nf*nt*params.nbits/8;
	uint8_t* in_buffer_device;
	printf("Copy in buffer size = %d MB per ant = %d MB TOTAL \n", in_buffer_bytes_per_ant/(1024l*1024l), in_buffer_bytes_per_ant*source->nants()/(1024l*1024l));
	gpuErrchk( cudaMalloc((void**) &in_buffer_device, in_buffer_bytes_per_ant*source->nants() ));


	DadaSink* dada_sink = NULL;
	if (params.export_dada_key != -1) {
		char* hdr = NULL;
		if (dada_source != NULL) {
			hdr = dada_source->get_source_at(0)->get_header();
		}
		dada_sink = new DadaSink(params.export_dada_key, hdr, params);
	}

	const bool dump_data = params.dump_data;
	const int nbeams_out = params.nbeams_out;
	// rescale output buffer
	array4d_t rescale_buf;
	rescale_buf.nw = nbeams_out;
	rescale_buf.nx = nf;
	rescale_buf.ny = 1;
	rescale_buf.nz = nt;
	array4d_malloc(&rescale_buf, dump_data, true);

	// rescale junk buffer for first integration only - bleah
	array4d_t rescale_junk_buf;
	rescale_junk_buf.nw = nbeams_out;
	rescale_junk_buf.nx = nf;
	rescale_junk_buf.ny = 1;
	rescale_junk_buf.nz = nt;
	array4d_malloc(&rescale_junk_buf, false, true);

	// FDMT output buffer
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
	rescale.target_stdev = 1.0/sqrt(params.nbeams_summed);
	rescale.decay_constant = 0.35 * params.decay_timescale / source->tsamp(); // This is how the_decimator.C does it, I think.
	rescale.mean_thresh = params.mean_thresh;
	rescale.std_thresh = params.std_thresh;
	rescale.kurt_thresh = params.kurt_thresh;
	rescale.flag_grow = params.flag_grow;
	rescale.dm0_thresh = params.dm0_thresh;
	rescale.cell_thresh = params.cell_thresh;
	rescale.invert_freq = (params.foff < 0);
	rescale.subtract_dm0 = params.subtract_dm0;
	rescale.nt = nt;
	rescale.nf = nf;
	rescale.nbeams_per_ant = params.nbeams_per_antenna;
	rescale.nants = source->nants();
	rescale.polsum = params.polsum;
	rescale.nbits = source->nbits();
	rescale.in_order = source->data_order();
	// set guess of initial scale and offset to dm0 thresholding works
	printf("Rescaling to mean=%f stdev=%f decay constant=%f mean/std/kurtosis/dm0/Cell thresholds: %0.1f/%0.1f/%0.1f/%0.1f/%0.1f grow flags by %d channels\n",
			rescale.target_mean,rescale.target_stdev,
			rescale.decay_constant,
			rescale.mean_thresh, rescale.std_thresh, rescale.kurt_thresh,
			rescale.dm0_thresh, rescale.cell_thresh,
			rescale.flag_grow);
	Rescaler* rescaler = new Rescaler(rescale);
	rescaler->set_scaleoffset(1.0f, 0.0f); // Just pass it straight through without rescaling

	float flag_freqs_mhz[] = {1111.0f, 1144.0f};
	int num_flag_freqs = sizeof(flag_freqs_mhz) / sizeof(float);
	for (int flagi = 0; flagi < num_flag_freqs; flagi++) {
		//float freq = source->fch1() + c * source->foff();
		float freq = flag_freqs_mhz[flagi];
		int channel = int(roundf((freq - source->fch1())/source->foff()));
		if (channel >= 0 && channel < nf) {
			printf("Flagging channel %d at frequency %f\n", channel, freq);
			rescaler->flag_channel(channel);
		}
	}

	// Create fdmt
	fdmt_t fdmt;
	printf("Creating FDMT fmin=%f fmax=%f nf=%d nd=%d nt=%d nbeams=%d nbeams_alloc=%d\n",
			params.fmin, params.fmax, nf, nd, nt, nbeams_out, params.nbeams_alloc);
	fdmt_create(&fdmt, params.fmin, params.fmax, nf, nd, nt, nbeams_out, params.nbeams_alloc, dump_data);

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

	CandidateList candidate_list(params.max_ncand_per_block);

	// measure bytes used
	size_t gpu_free_bytes, gpu_total_bytes;
	gpuErrchk(cudaMemGetInfo( &gpu_free_bytes, &gpu_total_bytes ));

	// add signal handler
	signal(SIGHUP, &handle_signal);
	signal(SIGINT, &handle_signal);
	signal(SIGTERM, &handle_signal);
	uint64_t num_flagged_beam_chans = 0;
	uint64_t num_flagged_times = 0;

	// Create streams - one for each antenan
	const int MAX_NANT = 72;
	cudaStream_t streams[MAX_NANT];
	assert(source->nants() <= MAX_NANT);
	for (int i = 0; i < source->nants(); i++) {
		gpuErrchk(cudaStreamCreate(&streams[i]));
		//streams[i] = 0;
	}

	while (true) {
		if (stopped) {
			printf("Stopped due to signal received\n");
			break;
		}
		if (blocknum >= params.max_nblocks) {
			break;
		}

		rescaler->reset(rescale_buf); // set output buffer to zero - each rescale update will add the result into the buffer

		fdmt.t_copy_in.start();

//#pragma omp parallel
		for(int iant = 0; iant < source->nants(); iant++) {
			// read samples from input - one antenna at a time.
			void* read_buf;
			int this_nt = source->read_samples_ant(&read_buf, iant);
			if (this_nt != nt) { // WE've run out of samples
				stopped = true;
				break;
			}
			// File is in TBF order
			// Output needs to be BFT order
			// Do transpose and cast to float on the way through using GPU
			uint8_t* this_ant_buffer = in_buffer_device + iant*in_buffer_bytes_per_ant;
			gpuErrchk(cudaMemcpyAsync(this_ant_buffer,
					read_buf, in_buffer_bytes_per_ant*sizeof(uint8_t), cudaMemcpyHostToDevice, streams[iant]));
			//tproc.start();
			//trescale.start();
			if (blocknum == 0 && params.num_rescale_blocks > 0) { // if first block rescale and update with no
				// flagging so we can work out roughly what the scales are
				// Send output to junk buffer - silly but will fix later
				// TODO: Remove junk buffer to save memory
				rescaler->update_and_transpose(rescale_junk_buf, this_ant_buffer, rescaler->noflag_options, iant, streams[iant]);

				// update scale and offset
				rescaler->update_scaleoffset(rescaler->noflag_options, iant, streams[iant]);
				if (params.do_dump_rescaler) {
					dump_rescaler(-1, rescaler);
				}
			}

			// this time we rescale with the flagging turned on
			rescaler->update_and_transpose(rescale_buf, this_ant_buffer, rescaler->options, iant, streams[iant]);
			//trescale.stop();
			//tproc.stop();
		}
		gpuErrchk(cudaDeviceSynchronize()); // Synchonize after doing all those asynchronous, multistream things
		fdmt.t_copy_in.stop();

		if (stopped) {// if we've run out of samples
			break;
		}

		if (dump_data) {
			dumparr("inbuf", iblock, &rescale_buf);
		}
		// Do asynchronous copy to dada output using the copy stream for antenna 0
		if (dada_sink != NULL) {
			void* outptr = dada_sink->open_block();
			gpuErrchk(cudaMemcpyAsync(outptr,
					rescale_buf.d_device,
					array4d_size(&rescale_buf)*sizeof(rescale_dtype),
					cudaMemcpyDeviceToHost,
					streams[0]));
		}

		// Count how many times were flagged
		assert(params.num_rescale_blocks >= 0);
		array4d_copy_to_host(&rescaler->nsamps); // must do this before updaing scaleoffset, which resets nsamps to zero
		tproc.start();

		for(int i = 0; i < nf*params.nbeams_in_total; ++i) {
			int nsamps = (int)rescaler->nsamps.d[i]; // nsamps is the number of unflagged samples from this block
			int nflagged = rescaler->sampnum - nsamps;
			// rescale.sampnum is the total number of samples that has gone into the rescaler since resetting
			assert (nflagged >= 0);
			num_flagged_times += nflagged;
		}

		// do rescaling if required
		if (params.num_rescale_blocks > 0 && blocknum % params.num_rescale_blocks == 0) {
			for(int iant = 0; iant < source->nants(); ++iant) {
				rescaler->update_scaleoffset(rescaler->options, iant);
			}

			// Count how many  channels have been flagged for this whole block
			// by looking at how many channels have scale==0
			array4d_copy_to_host(&rescaler->scale);
			for(int i = 0; i < nf*params.nbeams_in_total; ++i) {
				if (rescaler->scale.d[i] == 0) {
					// that channel will stay flagged for num_rescale_blocks
					num_flagged_beam_chans += params.num_rescale_blocks;
				}

				// it looks here like I'm counting twice, as we increment num_flagged_times outside the rescale_blocks_guard
				// but I'm a bit wary here, bcasue of teh danger, danger
//				// Count how many times have been flagged for this block
//				// TODO: DANGER DANGER! This doesn't count flagged times if num_rescale_blocks = 0
//				// This gave me a long headache at LAX when I set -s 1e30 stupidly.
//				int nsamps = (int)rescaler->nsamps.d[i];
//				// nsamps is the number of unflagged samples in nt*num_rescale_blocks samples
//				int nflagged = nt*num_rescale_blocks - nsamps;
//				assert (nflagged >= 0);
//				num_flagged_times += nflagged;
			}

			if (params.do_dump_rescaler) {
				dump_rescaler(iblock, rescaler);
			}
		}

		if (blocknum >= params.num_rescale_blocks) {
			/// Execute the FDMT
			fdmt_execute(&fdmt, rescale_buf.d_device, out_buf.d);
			if (dump_data) {
				dumparr("fdmt", iblock, &out_buf, false);
				dumparr("ostate", iblock, & fdmt.ostate, true);
			}
			//total_candidates += boxcar_threshonly(&out_buf, sampno, thresh, max_ncand_per_block, mindm, sink);
			tboxcar.start();
			boxcar_do_gpu (
					&fdmt.ostate,
					&boxcar_data,
					&boxcar_history,
					&boxcar_discards,
					params.thresh, params.max_ncand_per_block, params.mindm, params.maxbc, &candidate_list);
			tboxcar.stop();
			int ncand = candidate_list.copy_to_sink(sink);
			if (ncand >= params.max_ncand_per_block - 1) {
				num_candidate_overflow_blocks++;
			}
			total_candidates += ncand;
			if (dump_data) {
				dumparr("boxcar", iblock, &boxcar_data, true);
			}
		}
		tproc.stop();

		// release dada block from output -
		if (dada_sink != NULL) {
			gpuErrchk(cudaStreamSynchronize(streams[0]));
			dada_sink->close_block();
		}

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

	const int nbeams_in_total = params.nbeams_in_total;
	double boxcar_ngops = (double)nbeams_out*(double)nt*(double)nd*2.0*(double)NBOX/1e9;
	double data_nsecs = blocknum*nt*source->tsamp();

	double flagged_percent = ((double) num_flagged_beam_chans) / ((double) nf*nbeams_in_total*blocknum) * 100.0;
	double dm0_flagged_percent = ((double) num_flagged_times) / ((double) blocknum*nbeams_in_total*nt*nf) * 100.0;
	cout << " FREDDA Finished" << endl;
	cout << "Found " << total_candidates << " candidates" << endl;
	cout << "Discarded " << total_discards << " candidates for being too wide."<< endl;
	cout << num_candidate_overflow_blocks << " blocks overflowed the candidate buffer"<<endl;
	cout << "Processed " << blocknum << " blocks = "<< blocknum*nt << " samples = " << data_nsecs << " seconds" << " at " << data_nsecs/tall.wall_total()<< "x real time"<< endl;
	cout << "Freq auto-flagged " << num_flagged_beam_chans << "/" << (nf*nbeams_in_total*blocknum) << " channels = " << flagged_percent << "%" << endl;
	cout << "DM0 auto-flagged " << num_flagged_times << "/" << (blocknum*nbeams_in_total*nt*nf) << " samples = " << dm0_flagged_percent << "%" << endl;
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
	if(dada_sink != NULL) {
		delete dada_sink;
	}
}

