#!/usr/bin/env python
"""
Convvert packet dump data from a packtdump file to a filterbank file. See craft_pktdump.py

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'
import logging
import sys

import askap.time
import numpy as np
import time
import pylab
import pickle as pickle
import datetime
import pytz
import warnings
from askap.craft.sigproc import SigprocFile
from askap.craft.crafthdr import DadaHeader
from askap.craft.freqconfig import FreqConfig
import askap.craft.pktdump as pktdump

        
def bat2mjd(bat):
    utcdt = askap.time.bat2utcDt(int(bat))
    mjd = askap.time.utcDt2mjd(utcdt)
    print('BAT conversion', bat, utcdt, mjd)
    return mjd

def get_channel(hdr, address):
    fpga = int(np.log2(hdr['srcDevice']))
    card = int(address[0].split('.')[-1]) -1
    assert fpga >= 0 and fpga <= 5, 'Invalid FPGA %s' % fpga
    assert card >=0 and card <= 8, 'Invalid card %s' % card

    channel = card*6 + fpga
    return channel

class OutputReorderBuffer(object):
    def __init__(self, nslots, slotsize, fout, dtype=np.float32):
        assert nslots > 0
        self._nslots = nslots # number of slotts in the buffer

        assert fout is not None
        self._fout = fout # output file
        self._currslot = 0 # current slot number
        self._integ_nums = np.ones(nslots, dtype=np.int64)*-1         # Integration number for each slot. -1 means empty
        self._data = np.zeros((nslots, slotsize), dtype=dtype) # Buffer
        self._num_late = 0 # Counter for number of late packets

    def get_currint(self):
        return self._integ_nums[self._currslot]

    def new_slot(self, newint):
        currint = self.get_currint()
        newslot = (self._currslot + 1) % self._nslots

        logging.debug('Allocating new slot. currslot %s newslot %s currint %s newint %s', self._currslot, newslot, currint, newint)

        if currint > 0:
            assert newint == currint + 1, 'Non consecuitive integrations newint: %s curint: %s' % (newint, currint)

        if self._integ_nums[newslot] >= 0: # if the new slot is not empty
            logging.debug('Writing slot %s to fout', newslot)
            data_tofile = self._data[newslot, :]
            num_nan = np.sum(np.isnan(data_tofile))
#            if num_nan > 0:
#                logging.info('Slot %s for integ %s has %s NaNs (currint %s)', newslot, (currint - self._nslots+1), num_nan, currint)
            self._fout.write(data_tofile)
            
        self._data[newslot, :] = np.nan # clear
        
        # The new slot's integration is the old one + 1
        self._integ_nums[newslot] = newint

        # update current slot number
        self._currslot = newslot

        return newslot
            
            
    def set_currint(self, intno):
        currint = self.get_currint()

        if intno > currint: # WE have a new integration on our hands
            # Allocate new slots
            numnew = intno - currint
            if numnew > self._nslots:
                raise ValueError('Probably invalid intno >> nslots: intno=%s nslots=%s numnew=%s' % (intno, self._nslots, numnew))
            for i in range(numnew):
                self.new_slot(currint + i + 1)

        # Get currint again, as allocating slots will change currint
        return self.get_currint()


    def put(self, intno, idxs, intdata):

        currint = self.set_currint(intno)
        intdiff = currint - intno
        assert intdiff >= 0, 'Invalid intidff. currint {} intno {} intdiff {}'.format(currint, intno, intdiff)

        if intdiff >= self._nslots: # this packet has fallen off the back off the slots. SOrry.
            self._num_late += 1
            logging.info('Late packet. currint %s intno %s intdiff %s num_late %s', currint, intno, intdiff, self._num_late)
            return None
        
        slotno = (self._currslot - intdiff) % self._nslots
        logging.debug('Saving intdata intno %s to slot %s with intno %s ', intno, slotno, self._integ_nums[slotno])
        assert self._integ_nums[slotno] == intno, 'Logic error. integ_nums %s slotnot %s integ_nums[slotno] = %s intno %s' % (self._integ_nums, slotno, self._integ_nums[slotno], intno)

        self._copy_data(slotno, idxs, intdata)

    def _copy_data(self, slotno, idxs, intdata):
        self._data[slotno, idxs] = intdata

    def flush(self):
        currint = self.get_currint()
        for i in range(self._nslots):
            newint = currint + 1
            self.new_slot(newint)

class DataWriter(object):
    def __init__(self, main_hdr,  values):
        self.outfiles = []
        self.sigproc_files = []
        if values.source_name is None:
            source_name = main_hdr.get('SOURCE', 'UNKNOWN')
        else:
            source_name = values.source_name

        for i in range(values.nfiles):
            sigproc_hdr = {}
            sigproc_hdr['tstart'] = float(main_hdr['TSTART'])
            sigproc_hdr['tsamp'] = float(main_hdr['TSAMP'])
            sigproc_hdr['nifs'] = 1
            sigproc_hdr['nchans'] = int(main_hdr['NCHAN'])
            sigproc_hdr['source_name'] = source_name
            sigproc_hdr['telescope_id'] = 7
            sigproc_hdr['fch1'] = float(main_hdr['FREQ'])
            sigproc_hdr['foff'] = float(main_hdr['BW'])
            sigproc_hdr['rawdatafile'] = values.file
            sigproc_hdr['nbits'] = int(main_hdr['NBIT'])
            sigproc_hdr['machine_id'] = 1
            sigproc_hdr['data_type'] = 1 # fitlerbank
            sigproc_hdr['az_start'] = 0.
            sigproc_hdr['barycentric'] = 1
            sigproc_hdr['pulsarcentric'] = 0
            sigproc_hdr['nsamples'] = 1
            sigproc_hdr['za_start'] = 0.0
            sigproc_hdr['src_raj'] = 0.0 # TODO
            sigproc_hdr['src_dej'] = 0.0 # TODO
            fout = '{}_f{:02}.fil'.format(values.file, i)
            logging.debug('Creating filterbank file %s with header %s', fout, sigproc_hdr)
            fout = SigprocFile(fout, 'wr+', sigproc_hdr)
            self.outfiles.append(fout.fin)
            self.sigproc_files.append(fout)

        self.elements_per_file = int(main_hdr['NCHAN'])
        self.main_hdr = main_hdr
        self.values = values

    def write(self, data):
        n = self.elements_per_file
        for ifile, outfile in enumerate(self.outfiles):
            data_for_file = data[ifile*n:(ifile+1)*n]
            data_for_file.tofile(outfile)



class InputDefragmenter(object):
    def __init__(self, fragment_lengths):
        self._frag_lengths = fragment_lengths
        self._nfrag = len(self._frag_lengths)
        self._frags = {}
        self._seqnos = {}

    def put(self, hdr, addr, data):
        key = (addr, hdr['srcDevice'])
        fragments = self._frags.get(key, {})
        seqno = hdr['packetSequence']
        currseqno = self._seqnos.get(key, seqno)
        fragid = hdr['fragmentID']
        expected_fraglength = self._frag_lengths[fragid]
        #assert len(data) == expected_fraglength, 'Unexpected data length. Was: {} expected {} fragid {}'.format(len(data), expected_fraglength, fragid)
        if len(data) != expected_fraglength:
            warnings.warn('Unexpected data length. Was: {} expected {} fragid {}'.format(len(data), expected_fraglength, fragid))
            return None

        logging.debug('Defrag packet key %s nfrags %s, currseqno %s, seqno %s fragid %s hdr=%s', key, len(fragments), currseqno, seqno, fragid, hdr)

        if seqno != currseqno: 
            logging.info('Sequence number skipped key %s. Expected seqno %s actual seqno  %s fragid %s', key, currseqno, seqno, fragid)
            self._frags[key] = {}
            self._seqnos[key] = seqno
        
        
        fragments = self._frags.get(key, {})
        fragments[fragid] = (hdr, addr, data)
        
        if len(fragments) == self._nfrag:
            buf = ''
            for i in range(self._nfrag):
                d = fragments[i][2]
                buf += d

            new_seqno = (seqno + 1) % 256

            logging.debug('Got all fragments for key %s. They were: %s. New seqno %s', key, list(fragments.keys()), new_seqno)

            self._frags[key] = {}
            self._seqnos[key] = new_seqno

            return buf
        else:
            self._frags[key] = fragments
            logging.debug('Saved fragment for key %s. Now have %s fragments: %s', key, len(fragments), list(fragments.keys()))

            return None


def _main():
    from argparse import ArgumentParser
    parser = ArgumentParser(description='converst packet dump data to filterbank file')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s','--source-name', help='Set source name in filterbank', default='unknown')
    parser.add_argument('--print-times', action='store_true', help='Print packet time stamps')
    parser.add_argument('--freqconfig', help='Frequency mapping headerfile')
    parser.add_argument('-o', '--outfile', help='Ouptut file base')
    parser.add_argument('-n','--nfiles',help='Number of files to write', default=72, type=int)
    parser.add_argument(dest='file')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)


    fin = values.file
    f = pktdump.pktopen(fin)
    main_hdr = f.hdr

    print(main_hdr)
    intcount = int(main_hdr['INT_CYCLES'])
    intime = float(main_hdr['INT_TIME'])
    
    batint = intime*27./30.
    
    ntimes = intcount*4
    first_bat = None
    #for (hdr, address, data) in pickle.load(f):
    nrx = 0
    currbat = 0
    MAX_INT_COUNT = 8
    frame1 = None
    bat1 = None
    #ignore_cards = ['10.2.13.1']
    ignore_cards = []
    last_times = {}

    cards = list(map(int, main_hdr['CARDS'].split(',')))
    ncards = len(cards)
    nchans_per_fpga = 8 
    nfpgas_per_card = 6
    nchan = ncards*nchans_per_fpga*nfpgas_per_card
    nbeams = int(main_hdr['NBEAM'])
    npols =  int(main_hdr['NPOL'])
    nbeampols = nbeams*npols
    dtypes = {32: np.float32, 8: np.uint8}
    dtype = dtypes[int(main_hdr['NBIT'])]
    nslots = MAX_INT_COUNT * 4
    slotsize = nchan * nbeams * npols
    assert nslots >= intcount
    
    fraglengths = {}
    MAX_INTCOUNT = 7
    MAX_PKT_NBYTES = 8188
    
    for ic in range(1, MAX_INTCOUNT+1):
        # nbytes = 72 beams * 8 channels * 4 bytes/point + header(4 timestamps * 8 bytes)
        nbytes = (72*8*4 + 4*8) * ic
        npkts = nbytes/MAX_PKT_NBYTES
        fraglen = [MAX_PKT_NBYTES]*npkts + [nbytes % MAX_PKT_NBYTES]
        assert sum(fraglen) == nbytes
        fraglengths[ic] = fraglen
        print('Fragments for intcount', ic, fraglen)

    defrag = InputDefragmenter(fraglengths[intcount])

    if values.freqconfig is not None:
        hdr = DadaHeader.fromfile(values.freqconfig)
        freq_config = FreqConfig.load_from_dada_header(hdr, cards=cards)
        logging.info('Frequency config loaded: %s', freq_config)
        nchan = freq_config.nchan_span # override nchan
        bw = freq_config.bw
        freq = freq_config.freq
    elif hasattr(f, 'dada_hdr'):
        freq_config = FreqConfig.load_from_dada_header(f.dada_hdr)
        logging.info('Loading frequency config from DADA header %s', freq_config)
        nchan = freq_config.nchan_span # override nchan
        bw = freq_config.bw
        freq = freq_config.freq
    else:
        bw = main_hdr['BW']
        freq = main_hdr['FREQ']
        freq_config = None

    main_hdr['BW'] = str(bw)
    main_hdr['FREQ'] = str(freq)


    for hdr, address, pktdata, pktdate in f:
        if address[0] in ignore_cards:
            logging.info('Skipping bad card ' + address[0])
            continue

        if hdr['packetType'] != 129:
            continue

        data = defrag.put(hdr, address, pktdata)
        if data is None:
            logging.debug('Skipping fragment %s', hdr)
            continue

        times = np.fromstring(data[0:8*2*MAX_INT_COUNT], dtype=np.uint64)
        data  = np.fromstring(data[ntimes*8:], dtype=dtype)
        d = times
        nts = ntimes/2
        start_frames = d[0:nts:2]
        stop_frames = d[1:nts:2]
        start_bats = d[nts:2*nts:2]
        stop_bats = d[nts+1:2*nts:2]
        pktbat = start_bats[0]
        pktframe = start_frames[0]
        hdrbat = hdr['bat']
        dev = hdr['srcDevice']
        fpga = int(round(np.log2(dev)))
        seqno = hdr['packetSequence']
        bf_no = int(address[0].split('.')[3])
        ibf = cards.index(bf_no)

        if values.print_times:
            print('times', times)
            print('start bats', start_bats)
            print('stop bat', stop_bats)
            print('bat diffs', stop_bats - start_bats)
            print('start frames', start_frames)
            print('stop frames', stop_frames)
            print('frame diffs', stop_frames - start_frames)

        chan = get_channel(hdr, address)
        #print start_frames - stop_frames, type(start_frames), start_frames.dtype, start_frames, (start_frames - stop_frames), (start_frames + stop_frames), start_frames[0] + stop_frames[0]
        
        fdiff = stop_frames - start_frames 
        bdiff = stop_bats - start_bats 

        
        if first_bat is None:
            first_bat = pktbat + 5e4
            currbat = first_bat
            # This BAT doesn't seem to convert nicely to an mjd.
            print('Scheduled for first_bat', first_bat, pktbat)
            continue

        if pktbat < first_bat:
            print('Before first bat. Skipping ', first_bat, pktbat, first_bat - pktbat, hdr, len(data), d)
            continue

        if frame1 is None:
            frame1 = pktframe
            bat1 = pktbat
            hbat1 = hdrbat
            currint = None
            timezone = pytz.timezone('Australia/Perth')
            pktdate_awst = pktdate.replace(tzinfo=timezone)
            print('date', pktdate, pktdate_awst)
            #BAD IDEA! Isn't synced with data, and is for the wrong packet! But anway.
            # TODO: uess BAT as its missingthe top 8 (or is it 16???) bits
            main_hdr['TSTART'] = str(askap.time.utcDt2mjd(pktdate_awst)) 
            fout = DataWriter(main_hdr, values)
            outbuf = OutputReorderBuffer(nslots, slotsize, fout)

        if pktframe < frame1:
            continue

        # TODO: pktframes incremetn by < 1 integration over FPGAs. Grrr.
        # Remove INT when we're able
        pkt_frame_diff =  pktframe - frame1

        intno1 = int(np.floor(pkt_frame_diff/intime)) # assumes pktframe is right to within an integration 
        intno2 = int(np.floor(pkt_frame_diff/(intime*intcount)))*intcount # essumes pktframe is right to within an integration*incount
        # INTNO2 DOESN"T WORK!! See CRAFT-6

        # the actual intno used throughout the rest of the code - this is the integration number of the first integration in the defragmented packet
        intno = intno1

        if intno % 1000 == 0:
            print('INTNO', intno, pktframe, frame1, intime, pktframe-frame1, type(pktframe), type(frame1))

        intno_bat = (pktbat - bat1) / batint

        last_pktbat, last_pktframe, last_hdrbat, last_seqno, last_pktdate  = last_times.get((address, dev), (0, 0, 0, 0, None))

        datazero = np.all(data == 0)
        if last_pktdate is not None and pktdate is not None:
            pktdiff = str((pktdate - last_pktdate).microseconds/1e3)
            pktdatestr = str(pktdate.time())
        else:
            pktdiff = ''
            pktdatestr = ''

        if values.verbose:
            print('\t'.join(map(str, (address[0], hdr['srcDevice'],
#                                  pktdatestr, pktdiff,
                                  seqno, seqno-last_seqno, hdr['fragmentID'],
                                  hdrbat, hdrbat - hbat1, hdrbat - last_hdrbat,
                                  pktbat, pktbat-bat1,  pktbat-last_pktbat, '%0.3f' % intno_bat,
                                  pktframe, pktframe - frame1,pktframe - last_pktframe, '%0.3f' % intno,
                                  'fdiff', fdiff, bdiff, chan, datazero, len(data)
                                  ))))

        last_times[(address, dev)] = (pktbat, pktframe, hdrbat, seqno, pktdate)

        # for i10000s3 this makes the same channel live in the same time slot
        #data.shape = (intcount, nchans_per_fpga, nbeams)
        data.shape = (nchans_per_fpga, intcount, nbeampols)

        if freq_config is None:
            fidxs = np.arange(chan*nchans_per_fpga, (chan+1)*nchans_per_fpga)
        else:
            fidxs = freq_config.chanmaps[ibf, fpga*nchans_per_fpga:(fpga+1)*nchans_per_fpga]

        for i in range(intcount):
            # only writes 1 IF
            for b in range(nbeampols):
                d = data[:, i, b]
                idxs = fidxs + b*nchan
                outbuf.put(intno + i, idxs, d)

if __name__ == '__main__':
    _main()
