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
#include "ascii_header.h"

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
			"   -M M - Channel Mean relative change threshold. e.g. 1% is -M 0.01 \n"
			"   -T T - Channel StdDev relative changed flagging threshold. e.g. 30% is -T 0.3 \n"
			"   -K K - Channel Kurtosis threshold (3 is pretty good)\n"
			"   -G G - GTEST threshold - 0.25 is good. - Also must specify -I\n"
			"   -I I - Number of samples per integration - required if -G is specified\n"
			"   -z Z - Zap times with 0 DM above threshold Z\n"
			"   -C C - Zap time/frequency cells with S/N above threshold C\n"
			"   -F FILE - Flag frequencies contained in this file\n"
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

void FreddaParams::parse(int _argc, char* _argv[]) {
	argc = _argc;
	argv = &_argv[0];
	printf("Fredda version %s starting. Cmdline: ", VERSION);
	for (int c = 0; c < argc; ++c) {
		printf("%s ", argv[c]);
	}
	printf("\n");
	int ch;
	while ((ch = getopt(argc, argv, "d:t:s:o:x:r:S:B:DRg:M:T:U:K:G:I:C:F:n:m:b:z:N:X:uhp")) != -1) {
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
			//flag_grow = atoi(optarg);
			gtest_thresh = atof(optarg);
			break;
		case 'I':
			nsamps_per_int = atoi(optarg);
			break;
		case 'C':
			cell_thresh = atof(optarg);
			break;
		case 'F':
			flag_file = optarg;
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
	assert(gtest_thresh > 0);

	// Either INFINITY, -1 or non-infinity and defined are OK - this is exclusive or
	if ((gtest_thresh != INFINITY) && nsamps_per_int == -1) {
		printf("If you specify -G, you also need to specify -I\n");
		exit(EXIT_FAILURE);
	}
}

void FreddaParams::set_source(DataSource& _source) {
	// rescaling always puts the minimum frequency as the channel 0
	source = &_source;
	if (source->foff() > 0) {
		out_freq = source->fch1();
	} else {
		out_freq = source->fch1() + source->nchans()*source->foff();
	}
	assert(out_freq <= source->fch1());
	out_foff = abs(source->foff());

	nbeams_per_antenna = source->nbeams()*source->npols(); // number of beams including polarisations
	nbeams_in_total = nbeams_per_antenna*source->nants();
	npols_in = source->npols();
	if (polsum) { // assume polsum and antsum
		nbeams_out = source->nbeams();
		npols_out = 1;
		assert(nbeams_per_antenna %2 == 0);
	} else { // ant sum only
		nbeams_out = source->nbeams()*source->npols();
		npols_out = source->npols();
	}
	nbeams_summed = (float(nbeams_in_total)/float(nbeams_out));
	nf = source->nchans();
	nbits = source->nbits();
	foff =  (float) source->foff();
	fmax = (float) source->fch1() - foff; // The FDMT seems to want this offset to make sense of the world. Not sure why.
	fmin = fmax + nf*foff;

	if (nd < 0) { // Flip the band to calculate negative DMs
		nd = -nd; // make nd positive -otherwise array sizes get confuddled
		// FDMT requres fmin < fmax
		// rescaling will invert the channels now that we've changed the sign of foff
		foff = -foff;
	}
}

void FreddaParams::to_dada(char header_buf[])
{

	ascii_header_set(header_buf, "NT", "%d", nt);
	ascii_header_set(header_buf, "ND", "%d", nd);
	ascii_header_set(header_buf, "SEEK_SECONDS", "%f", seek_seconds);
	ascii_header_set(header_buf, "NUM_RESCALE_BLOCKS", "%d", num_rescale_blocks);
	ascii_header_set(header_buf, "DECAY_TIMESCALE", "%f", decay_timescale);
	ascii_header_set(header_buf, "THRESH", "%f", thresh);
	ascii_header_set(header_buf, "OUT_FILENAME", "%s", out_filename);
	ascii_header_set(header_buf, "DUMP_DATA", "%d", dump_data);
	ascii_header_set(header_buf, "DO_DUMP_RESCALER", "%d", do_dump_rescaler);
	ascii_header_set(header_buf, "CUDA_DEVICE", "%d", cuda_device);
	ascii_header_set(header_buf, "KURT_THRESH", "%f", kurt_thresh);
	ascii_header_set(header_buf, "STD_THRESH", "%f", std_thresh);
	ascii_header_set(header_buf, "MEAN_THRESH", "%f", mean_thresh);
	ascii_header_set(header_buf, "DM0_THRESH", "%f", dm0_thresh);
	ascii_header_set(header_buf, "CELL_THRESH", "%f", cell_thresh);
	ascii_header_set(header_buf, "FLAG_GROW", "%d", flag_grow);
	ascii_header_set(header_buf, "MAX_NCAND_PER_BLOCK", "%d", max_ncand_per_block);
	ascii_header_set(header_buf, "MINDM", "%d", mindm);
	ascii_header_set(header_buf, "MAXBC", "%d", maxbc);
	ascii_header_set(header_buf, "MAX_NBLOCKS", "%d", max_nblocks);
	ascii_header_set(header_buf, "NBEAMS_ALLOC", "%d", nbeams_alloc);
	ascii_header_set(header_buf, "SUBTRACT_DM0", "%d", subtract_dm0);
	ascii_header_set(header_buf, "POLSUM", "%d", polsum);
	ascii_header_set(header_buf, "UDP_HOST", "%s", udp_host);
	ascii_header_set(header_buf, "UDP_PORT", "%d", udp_port);
	ascii_header_set(header_buf, "EXPORT_DADA_KEY", "%x", export_dada_key);
	ascii_header_set(header_buf, "OUT_FREQ", "%f", out_freq);
	ascii_header_set(header_buf, "OUT_FOFF", "%f", out_foff);
	ascii_header_set(header_buf, "NBEAMS_PER_ANTENNA", "%d", nbeams_per_antenna);
	ascii_header_set(header_buf, "NBEAMS_IN_TOTAL", "%d", nbeams_in_total);
	ascii_header_set(header_buf, "NPOLS_IN", "%d", npols_in);
	ascii_header_set(header_buf, "NBEAMS_OUT", "%d", nbeams_out);
	ascii_header_set(header_buf, "NPOLS_OUT", "%d", npols_out);
	ascii_header_set(header_buf, "NBEAMS_SUMMED", "%f", nbeams_summed);
	ascii_header_set(header_buf, "NF", "%d", nf);
	ascii_header_set(header_buf, "NBITS", "%d", nbits);
	ascii_header_set(header_buf, "FOFF", "%f", foff);
	ascii_header_set(header_buf, "FMAX", "%f", fmax);
	ascii_header_set(header_buf, "FMIN", "%f", fmin);
	ascii_header_set(header_buf, "SOURCE_FCH1", "%f", source->fch1());
	ascii_header_set(header_buf, "SOURCE_FOFF", "%f", source->foff());
	ascii_header_set(header_buf, "SOURCE_TSTART", "%0.12f", source->tstart());
	ascii_header_set(header_buf, "SOURCE_TSAMP", "%0.12f", source->tsamp());
	ascii_header_set(header_buf, "SOURCE_NANTS", "%d", source->nants());
	ascii_header_set(header_buf, "SOURCE_NBEAMS", "%d", source->nbeams());
	ascii_header_set(header_buf, "SOURCE_NPOLS", "%d", source->npols());
	ascii_header_set(header_buf, "SOURCE_NCHANS", "%d", source->nchans());
	ascii_header_set(header_buf, "SOURCE_NBITS", "%d", source->nbits());
	ascii_header_set(header_buf, "SOURCE_NSAMPS_PER_INT", "%d", source->nsamps_per_int());
	ascii_header_set(header_buf, "SOURCE_DATA_ORDER", "%d", source->data_order());
	ascii_header_set(header_buf, "SOURCE_NAME", "%s", source->name());
	ascii_header_set(header_buf, "SOURCE_ANTENNA_NAME", "%s", source->antenna_name());
}

