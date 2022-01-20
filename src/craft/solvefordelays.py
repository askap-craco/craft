#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2015
"""
import pylab
import glob
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
from .crafthdr import DadaHeader

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-s','--show', action='store_true', default=False, help='Show plots')
    parser.add_argument('-t','--scan', help='Scan to use, 1st=0, 2nd=1, ...', type=int, default=0)
    parser.add_argument('-p','--hfiles', help='Path to header files', type=str, required=True, nargs='*')
    parser.add_argument('-d','--data', help='Path to .real and .imag data files', dest='files', required=True, nargs='*')
    parser.add_argument('-b','--beam', help='Beam number', type=int, required=True)
    parser.add_argument('-i','--sbid', help='SB ID', type=int, required=True)
    values = parser.parse_args()

    beamno	  = values.beam
    mode	  = '?'
    sb		  = 'SB' + str(values.sbid)
    scanid	  = values.scan
    downlbeam	  = ''
    fileextension = 'vcraft.hdr'

    # get antenna range from the distinct directories in the path, 
    # and assume:
    # - cards range = determine from the available header files
    # - fpga range = 0-5
    # - nfpgachans = 8 

    nfpgas	= 6
    nfpgachans	= 8

    SN_THRESHOLD = 2.0

    #hfiles = glob.glob('{}ak*/beam{}/ak*_c*_f*.{}'.format(hpath, beamno, fileextension))
    hfiles = values.hfiles
    hdrs = [DadaHeader.fromfile(f) for f in hfiles]
        
    antennalist = sorted(set([f['ANT'][0].lower() for f in hdrs]))
    print('Antennalist: ', antennalist)

    # determine the available cards from the header files
    cardsallhfiles = []
    for f in hfiles:
        hdr = DadaHeader.fromfile(f)
    	card = int(hdr.get_value('CARD_NO'))
	cardsallhfiles.append(card)
    cardnumbers = np.unique(cardsallhfiles)
    print('Cardnumbers: ', cardnumbers)

    ncards = len(cardnumbers)


    # fill hdrdata array with data from header files
    hdrdata = np.zeros((len(antennalist), ncards, nfpgas, nfpgachans), dtype=int)
    for f in hfiles:
        hdr = DadaHeader.fromfile(f)

	antenna_s = 'ak{:02d}'.format(int(hdr.get_value('ANTENNA_NO')))
    	antenna_i = antennalist.index(antenna_s)
    
    	card = int(hdr.get_value('CARD_NO'))
	cardindex =  np.where(cardnumbers == card)
    	fpga = int(hdr.get_value('FPGA_ID'))
    	hdrdata[antenna_i, cardindex, fpga, :] = hdr.get_value('FREQS').split(',')

	mode = int(hdr.get_value('MODE'))
	downlbeam = int(hdr.get_value('BEAM'))
    	if downlbeam != beamno:
    	    print('Beam in header file {} ({}) does not match beam argument {}'.format(f, downlbeam, beamno))


    # check that each card-fpga combination handles the same frequencies on each antenna
    for f in range(nfpgas):
        for c in range(ncards):
            for a in range(1, len(antennalist)):
                if hdrdata[a-1, c, f, :].all() == hdrdata[a, c, f, :].all():
                    #print 'freqs for C{} F{} in antennas {} and {} are the same'.format(c, f, a-1, a)
                    pass
                else:
                    print('!!! freqs for C{} F{} in antennas {} and {} are NOT the same'.format(c, f, a-1, a))

    # make a mapping for frequency - card/fpga combination;
    # sorted from low freq to high freq
    tmpmap = []
    # use the reference antenna for the mapping
    antenna = 0
    for f in range(nfpgas):
        for c in range(ncards):
            for ch in range(nfpgachans):
                freq = hdrdata[antenna, c, f, ch]
                tmpmap.append([freq, c, f])

    hwmap = sorted(tmpmap, key=lambda t: t[0])


    nant = len(antennalist) - 1
    nchan = 54
    freqs = (ncards * nfpgas * nfpgachans) - 1

    dfiles = values.files
    if not dfiles:
    	print('No data files in {} for beam{}'.format(values.files, beamno))
	sys.exit()

    cf_sorted_spec = np.zeros((nant, len(dfiles), ncards * nfpgas * nfpgachans, nchan), dtype=np.complex64)
    cf_sorted_av   = np.zeros((nant, len(dfiles), ncards * nfpgas, nchan), dtype=np.complex64)
    timestamps = dfiles

    for f in dfiles:
    	scan = dfiles.index(f)
        rfile = f
        imagfile = f.replace('.real','.imag')
        rdata = np.loadtxt(rfile)
        idata = np.loadtxt(imagfile)
        assert rdata.shape == idata.shape
        xaxis_all = rdata[:, 0]
        xaxis_all.shape = (nant, -1) # should be the same for all antennas
        for a in range(nant):
            xaxis = xaxis_all[a, :]
            assert np.all(xaxis == xaxis_all[0, :])


        assert np.all(rdata[:, 0] == idata[:, 0])
        fulld = np.zeros(len(rdata), dtype=np.complex64)
        fulld.real = rdata[:, 1]
        fulld.imag = idata[:, 1]
        fulld.shape = (nant, -1, nchan)
        delayspec = np.fft.fftshift(np.fft.fft(fulld, axis=2), axes=2)
        
	# sort the data as a function of card-fpga
	for a in range(nant):
	    for i in range(delayspec.shape[1]):
	    	# frequencies are sorted from high to low in delayspec, so look up the
            	# frequency index in hwmap at freqs-i
	    	freqoffset = int(np.argwhere(hdrdata[a,hwmap[freqs-i][1],hwmap[freqs-i][2],:] == hwmap[freqs-i][0]))
	    	cardoffset = hwmap[freqs-i][1] * nfpgas * nfpgachans
	    	fpgaoffset = hwmap[freqs-i][2] * nfpgachans
	    	# bins are sorted as a function of card-fpga-frequency
	    	bin = cardoffset+fpgaoffset+freqoffset
            	cf_sorted_spec[a, scan, bin] = delayspec[a, i, :]


	# combine all 8 channels of each fpga
	for a in range(nant):

	    fpgabin = 0
	    for i in range(0, cf_sorted_spec.shape[2], nfpgachans):
		for j in range(nfpgachans):
		    cf_sorted_av[a, scan, fpgabin] += abs(cf_sorted_spec[a, scan,i+j,:].T)

		fpgabin += 1


    # create and fill array with delays for each fpga
    delays = np.zeros((len(antennalist), len(timestamps), ncards * nfpgas), dtype=int)
    for antenna in range(cf_sorted_av.shape[0]):
	for scan in range(cf_sorted_av.shape[1]):
	    for fpga in range(cf_sorted_av.shape[2]):
		offset = 0
	    	std = np.std(cf_sorted_av[antenna, scan, fpga, :])
		zerodelay = cf_sorted_av[antenna, scan, fpga, cf_sorted_av.shape[3]/2]
		mean = np.mean(cf_sorted_av[antenna, scan, fpga, :])
		sn = zerodelay / mean

		# check if the bin at zerodelay is significant
		if sn < SN_THRESHOLD:
		    offsetmax = np.amax(cf_sorted_av[antenna, scan, fpga, :])
		    snoffset = offsetmax / mean
		    ioffsetmax = np.argmax(cf_sorted_av[antenna, scan, fpga, :])
		    # check if the maximum is significant, otherwise leave offset = 0
		    if snoffset > SN_THRESHOLD:
		    	offset = cf_sorted_av.shape[3]/2 - ioffsetmax

                print('DELAYS', delays.shape)
		delays[antenna+1, scan, fpga] = offset
		#print '{} c{} f{} {}'.format(antennalist[antenna+1], fpga/nfpgas + 1, fpga%nfpgas, offset)

    # check whether the delays stay the same for each scan
    for antenna in range(delays.shape[0]):
    	for scan in range(1, delays.shape[1]):
	    if not np.array_equal(delays[antenna, 0, :], delays[antenna, scan, :]):
	    	print('ANTENNA {}: delays for scan {} differ with respect to first scan {}'.format(antennalist[antenna], timestamps[scan], timestamps[0]))

    if scanid > len(timestamps) - 1:
    	print('WARNING: Invalid scanid {}, there are only {} scans, using first scan.'.format(scanid, len(timestamps)))
    	scanid = 0

    # write the delays from the first scan to file
    outfilename = 'hwdelays_' + sb + '_' + str(timestamps[scanid]) + '.txt' 
    outfile = open(outfilename, 'w')
    for antenna in range(delays.shape[0]):
	for fpga in range(delays.shape[2]):
	    # get card number
	    card = cardnumbers[fpga/nfpgas]
	    outfile.write('{}_c{}_f{}.{} {}\n'.format(antennalist[antenna], card, fpga%nfpgas, fileextension, delays[antenna, scanid, fpga]))

    outfile.close()

    if values.show:
    	lowedge_cards = []
    	coffset = 3
    	for c in range(ncards):
	    lowedge_cards.append((c*nfpgas) + coffset)

    	# nrows = number of baselines per figure
	nrows = 3

	# determine number of figures
	if nant%nrows == 0:
	    nfigs = nant/nrows
	else:
	    nfigs = nant/nrows + 1

	for nfig in range(nfigs):
	    # if number of baselines is not a multiple of nrows and if this is the last figure,
	    # determine the number of rows
	    if nant%nrows != 0 and nfig == nfigs - 1:
	    	rowsinfig = nant - (nfig*nrows)
	    else:
	    	rowsinfig = nrows
 
    	    fig, ax = pylab.subplots(rowsinfig, cf_sorted_av.shape[1])
    	    ax = ax.flatten()
    	    fig.subplots_adjust(left=0.03, bottom=0.01, right=1, top=0.93, wspace=0.05, hspace=0.1)

    	    extent = (0, ncards*nfpgas, -nchan/2, nchan/2)
    	    for a in range(rowsinfig):
	    	for t in range(cf_sorted_av.shape[1]):
	    	    # set yaxis label for left column, remove labels for all other timestamps
	    	    if t == 0:
		    	ax[a*cf_sorted_av.shape[1]].set_ylabel('delay (samples)')	    	
	    	    else:
    	            	ax[t+a*cf_sorted_av.shape[1]].yaxis.set_major_formatter(plt.NullFormatter())

	    	    # set title for top row
	    	    if a == 0:
	    	    	ax[t].set_title(timestamps[t])

	    	    # set xlabels
	            ax[t+a*cf_sorted_av.shape[1]].set_xticks(lowedge_cards, minor=False)
	            ax[t+a*cf_sorted_av.shape[1]].set_xticklabels(cardnumbers, minor=False)

		    # indicate cards with verical lines
	    	    for line in lowedge_cards:
	    	    	ax[t+a*cf_sorted_av.shape[1]].axvline(line - coffset, color='black')
            	    ax[t+a*cf_sorted_av.shape[1]].imshow(np.log10(abs(cf_sorted_av[(nfig*nrows)+a, t, :, :].T)),  interpolation='nearest', aspect='auto', extent=extent)

	    stitle = '{} / mode{} / beam{} / '.format(sb, mode, beamno)
	    for r in range(rowsinfig):
	        stitle += 'Row{}: {} - {} '.format(r+1, antennalist[0], antennalist[(nfig*nrows)+r+1])
    	    fig.suptitle('{}'.format(stitle))

        pylab.show()


if __name__ == '__main__':
    _main()
