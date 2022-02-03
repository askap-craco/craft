#!/usr/bin/env python
"""
Sets up craft and downlaods channel mapping

Copyright (C) CSIRO 2016
"""
__author__ = 'Keith Bannister <keith.bannister@csiro.au>'
import logging
import sys
import askap.craft.beamformer 
from askap.craft.beamformer import CraftBeamformer, FilterTypes, ZoomModes
from askap.craft.freqconfig import FreqConfig
import numpy as np
from multiprocessing import Pool

fpgaAddress = [0x1, 0x2, 0x4, 0x8, 0x10, 0x20, 0x40]
MAX_BF_CARD_BW = 48                      #@brief Processed BW of each beamformer card in non-zoom mode
MAX_BF_FPGA_BW = MAX_BF_CARD_BW/6        #@brief Processed BW of each beamformer fpga in non-zoom mode

def getCraftSubbands(bf):
    #m_craftSubbands[i] = (m_bullant->getFPGASubbandID(fpgaAddress[i]) & 0x7f) * MAX_BF_FPGA_BW;
    subbands = []
    for i in range(6):
        b = bf.getBullantCard().getFPGASubbandID(fpgaAddress[i] & 0x7f) * MAX_BF_FPGA_BW
        subbands.append(b)

    print('Got subbands', subbands)

    return subbands
        

'''
  m_craftSkyFreqs.resize(0);
  for (unsigned int id = 0; id < 6; id++) {
    for(unsigned int freq = 0; freq < 8; freq++) {
      float skyfreq = m_skyFreq[m_craftSubbands[id]+freq].freq;
      m_craftSkyFreqs.push_back(skyfreq);
    }
  }
'''
def getCraftSkyFrequencies2(bf):
    status = bf.getStatus()
    allfreqs =  [chan.freq for chan in status.channelNumbers]
    subbands = getCraftSubbands(bf)
    skyfreqs = []

    for _id in range(6):
        for _freq in range(8):
            skyfreq = allfreqs[subbands[_id] + _freq]
            skyfreqs.append(skyfreq)

    return skyfreqs

def _main():
    from argparse import ArgumentParser
    from askap.craft.cmdline import strrange
    parser = ArgumentParser(description='Sets up craft and downloads channel mapping')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-a', '--antennas', help='Antenna number, comma sparated or dash. e.g. "3", "13-15", "1,2,14"', type=strrange)
    parser.add_argument('-c','--cards', help='Card number. comma separated, or dash. e.g. "3", "1-7" or "1,2,4-7"', default="1-7", type=strrange)
    parser.add_argument('-i','--int-time', help='Integration time (1-65535)', default=1000, type=int)
    parser.add_argument('-s','--int-cycles', help='Number of cycles to combine (1-7)', default=7, type=int)
    parser.add_argument('-n','--beam-number', help='Beam number to save voltages for (0-35). Unspecified for all beams', default=None)
    parser.add_argument('-m','--buffer-mode', help='Buffer number of bits', choices=[16, 8, 4, 1], type=int, default=16)
    parser.add_argument('-b','--band',help='Band number 0 = 1200 MHz BPF, 1 = 1450 MHz BPF, 2 - 1800 MHz BPF, 3 - LPF', choices=[0,1,2,3], type=int)
    parser.add_argument('-z','--zoom',help='Zoom mode', choices=[0,1,2,3,4,5], default=0, type=int)
    parser.add_argument('-f','--center-freq', help='Center frequency (MHz)', type=int)
    parser.add_argument('-d','--pushd-delay', help='PUsh download delay -1=auto, 0=none, x=value (seconds)', type=float, default=0)
    parser.add_argument('--program-fpga', help="Force programming of FPGA. 0=Don't. 1=Do", type=int, choices=[0,1], default=0)
    parser.add_argument('-o','--freqfile', help='File to write frequency maping to', default='.freqs')
    parser.add_argument('--corrblock', help='Correlator block to which to apply center freq, unspecified or 0 for default', type=int, choices=[0,1,2,3,4,5,6,7,8], default=0)
    parser.add_argument('--num-threads', help='Number of threads to spawn', type=int, default=1)
    
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    pool = Pool(processes=values.num_threads)
    antv = [(ant, values) for ant in values.antennas]
    pool.map(setup_craft_ant, antv)


def setup_craft_ant(antv):
    ant, values = antv
    iant = values.antennas.index(ant)
    nant = len(values.antennas)

    int_time = int(values.int_time)
    assert int_time >=1 and int_time<= 65535, 'Invalid integration time {}'.format(int_time)
    
    int_cycles = int(values.int_cycles)
    assert int_cycles >=1 and int_cycles <= 7, 'Invalid int cycles {}'.format(int_cycles)

    bufmode_map = {16:0, 8:1, 4:2, 1:3}
    assert values.buffer_mode in list(bufmode_map.keys()), 'Invalid buffer mode'
    bufmode = bufmode_map[values.buffer_mode]

    if values.beam_number is None:
        beam_number = -1
    else:
        beam_number = int(values.beam_number)
        assert beam_number >= 0 and beam_number < 36, 'Invalid beam number to save {}'.format(beam_number)
        bufmode += 4 # This tells the firmware to save only 1 beam

    programFPGA = values.program_fpga  # please dont' program FPGA - just be nice
    bandNo = values.band
    centerFreq = values.center_freq
    zoomMode = values.zoom
    corrblock = values.corrblock

    bfs = []
    pushDownloadDelay = values.push_delay
    print('Looking at cards', values.cards)

    foutname  = 'ak{}_{}'.format(ant, values.freqfile)
    logging.info('Writing file %s', foutname)
    fout = open(foutname, 'w')
    fout.write('HDR_VERSION 1.0\n')
    fout.write('HDR_SIZE {}\n'.format(4096*16))
    fout.write('HDR DADA\n')
    fout.write('ANTNO {}\n'.format(ant))
    fout.write('PROGRAMFPGA {}\n'.format(programFPGA))
    fout.write('CORRBLOCK {}\n'.format(corrblock))
    fout.write('CENTERFREQ {}\n'.format(centerFreq))
    fout.write('ZOOMMODE {}\n'.format(zoomMode))
    fout.write('INT_CYCLES {}\n'.format(int_cycles))
    fout.write('INT_TIME {}\n'.format(int_time))
    fout.write('CRAFT_SETUP_BEAM_NUMBER {}\n'.format(beam_number))
    fout.write('CRAFT_SETUP_BUFMODE {}\n'.format(bufmode))
    fout.write('CRAFT_SETUP_PUSHDELAY {}\n'.format(pushDownloadDelay))

    
    fout.write('NUM_BEAMFORMERS {}\n'.format(len(values.cards)))
    
    fstart = 768
    packet_interval_sec = float(int_cycles*int_time)*32./27./1e6
    
    for icard, card in enumerate(values.cards):
        bf = CraftBeamformer(ant, card)
        bf.connect()
        total_card_no = icard + iant*nant
        
        
            #bf.check_craft_enabled()
            #bf.startupBeamformer(programFPGA, FilterTypes.values[bandNo], int(centerFreq), ZoomModes.values[zoomMode], pushDownloadDelay)
            #bf.getBullantCard().enableOffboardLinks()
        delay_sec = packet_interval_sec*float(total_card_no)
        delay_samps = int(np.ceil(delay_sec/8e-9))

        print('DELAYS', ant, iant, card, icard, packet_interval_sec, delay_sec, delay_samps)

        if values.push_delay == -1:
            bf.getBullantCard().setDownloadDelay(delay_samps)

        bf.setCenterFrequency(FilterTypes.values[bandNo], ZoomModes.values[zoomMode], centerFreq, programFPGA)
            #bf.setupCraft(int_time, int_cycles, bufmode, beam_number, 0)
        bf.wait_for_done()
        craftskyfreqs = bf.getCraftSkyFrequencies()
            #freqs = getCraftSkyFrequencies2(bf)
        f1 = fstart + 48*icard
        fakefreqs = np.arange(f1, f1+48)
        freqs = fakefreqs
        assert(len(freqs) == 48)
        
        status = bf.getStatus()
        allfreqs =  [chan.freq for chan in status.channelNumbers]
        myfreqs = allfreqs[48*icard:48*(icard+1)]
        finesteps =  [chan.fineStep for chan in status.channelNumbers]
        zooms =  [chan.zoom for chan in status.channelNumbers]
        
        fout.write('BEAMFORMER{}_CARDNO {}\n'.format(icard, card))
        fout.write('BEAMFORMER{}_ADDR {}\n'.format(icard, str(bf.data_addr)))
        fout.write('BEAMFORMER{}_CRAFTSKYFREQS {}\n'.format(icard, ','.join(map(str, craftskyfreqs))))
        fout.write('BEAMFORMER{}_MYFREQS {}\n'.format(icard, ','.join(map(str, myfreqs))))
        fout.write('BEAMFORMER{}_ALLSKYFREQS {}\n'.format(icard, ','.join(map(str, allfreqs))))
        fout.write('BEAMFORMER{}_FREQS {}\n'.format(icard, ','.join(map(str, freqs))))

            

if __name__ == '__main__':
    _main()
