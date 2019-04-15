/*
 * FreddaParams.cpp
 *
 *  Created on: 15 Apr 2019
 *      Author: ban115
 */

#include "FreddaParams.h"

#include <stdio.h>
#include <unistd.h>
#include <assert.h>

void runtest_usage() {
	fprintf(stderr,
			"cudafdmt [options] [infile [infile[ ...]]\n"
			"   -d D - Number of dispersion trials. Negative D computes negative DMs\n"
			"   -t T - Samples per block\n"
			"   -s S - Decay timescale\n"
			"   -o FILE - Candidate filename\n"
			"   -U host:port - UDP host:port to send candidates to\n"
			"   -x SN - threshold S/N\n"
			"   -D dump intermediate data to disk (SLOW)\n"
			"   -R dump rescaler data to disk\n"
			"   -B b - Process b beams simultaneously to save memory\n"
			"   -r R - Blocks per rescale update (0 for no rescaling)\n"
			"   -S S - Seek to this number of seconds before starting\n"
			"   -M M - Channel Mean relative change threshold (0.2 is OK)\n"
			"   -T T - Channel StdDev relative changed flagging threshold (0.2 is OK)\n"
			"   -K K - Channel Kurtosis threshold (3 is pretty good)\n"
//			"   -G N - Channel flag channel growing (flags N channels either side of a bad channel)\n"
			"   -z Z - Zap times with 0 DM above threshold Z\n"
			"   -C C - Zap time/frequency cells with S/N above threshold C\n"
			"   -u   - Subtract DM0 time series from spectrum\n"
			"   -p   - Sum polarisations\n"
			"   -n ncand - Maximum mumber of candidates to write per block\n"
			"   -m mindm - Minimum DM to report candidates for (to ignore 0 DM junk)\n"
			"   -b maxbc - Maximum boxcar to create a candidate. Candidates with peaks above this boxcar are ignored\n"
			"   -g G - CUDA device\n"
			"   -N N - Maximum number of blocks to process before quitting\n"
			"   -X x - Export incoherent sum data to this DADA key\n"
			"   -h Print this message\n"
			"    Version: %s\n"
			, VERSION);
	exit(EXIT_FAILURE);
}

FreddaParams::FreddaParams() {
	bzero(udp_host, 128);
}

FreddaParams::~FreddaParams() {
}

void FreddaParams::parse(int argc, char* argv[])
{
	printf("Fredda version %s starting. Cmdline: ", VERSION);
	for (int c = 0; c < argc; ++c) {
		printf("%s ", argv[c]);
	}
	printf("\n");
	int ch;
	while ((ch = getopt(argc, argv, "d:t:s:o:x:r:S:B:DRg:M:T:U:K:G:C:n:m:b:z:N:X:uhp")) != -1) {
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
		case 'R':
			do_dump_rescaler = true;
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
		case 'X':
			sscanf(optarg, "%x",&export_dada_key);
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
	assert(seek_seconds >= 0);
}

void FreddaParams::set_source(DataSource& source) {
	// rescaling always puts the minimum frequency as the channel 0
	if (source.foff() > 0) {
		out_freq = source.fch1();
	} else {
		out_freq = source.fch1() + source.nchans()*source.foff();
	}
	assert(out_freq <= source.fch1());
	out_foff = abs(source.foff());

	nbeams_per_antenna = source.nbeams()*source.npols(); // number of beams including polarisations
	nbeams_in_total = nbeams_per_antenna*source.nants();
	npols_in = source.npols();
	if (polsum) { // assume polsum and antsum
		nbeams_out = source.nbeams();
		npols_out = 1;
		assert(nbeams_per_antenna %2 == 0);
	} else { // ant sum only
		nbeams_out = source.nbeams()*source.npols();
		npols_out = source.npols();
	}
	nbeams_summed = (float(nbeams_in_total)/float(nbeams_out));
	nf = source.nchans();
	nbits = source.nbits();
	foff =  (float) source.foff();
	fmax = (float) source.fch1() - foff; // The FDMT seems to want this offset to make sense of the world. Not sure why.
	fmin = fmax + nf*foff;

	if (nd < 0) { // Flip the band to calculate negative DMs
		nd = -nd; // make nd positive -otherwise array sizes get confuddled
		// FDMT requres fmin < fmax
		// rescaling will invert the channels now that we've changed the sign of foff
		foff = -foff;
	}
}

