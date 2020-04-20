#!/usr/bin/env python
# coding: utf-8

"""
obtain antenna delays and phase
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import io
import os
import sys

import glob



class aipscor(object):
    def __init__(self, fringfile, scfile, bpfile, vcraft_dr=None, imfile=None, input_dr=None):
        self.vcraft_dr = vcraft_dr
        self.imfile = imfile
        self.input_dr = input_dr
        self.fringfile = fringfile #delays.sn.txt'
        self.scfile = scfile #selfcal.sn.txt'
        self.bpfile = bpfile #bandpass.bp.txt'
        if vcraft_dr is not None and imfile is not None:
            self.get_basic_info()
        
        
    def get_basic_info(self):
        # self.nant : total number of antennas
        # self.npol : total number of polarizations (should be either 1 or 2)
        # self.pol2beam : map of polarization to beam numbers
        # self.pols : list of given polarizations
        # self.beams : list of given beams
        # self.vfiles : list of vcraft files sorted by polarization
        # self.an_name : list of antenna names in akNN format
        # self.freqs : dictionary of frequencies for each FPGAs
        # self.nfreq : number of total frequency channels
        
        # number of antennas
        antennas = glob.glob(self.vcraft_dr+'/*/')
        self.nant = len(antennas)
        
        # number of polarizations
        pols = glob.glob(antennas[0]+'/*')
        self.npol = len(pols)
        
        # map polarizations to beam numbers
        self.pol2beam = {}
        self.pols = []
        self.beams = []
        for p in pols:
            self.pols += [p]
            beam = int(p.split('beam')[-1])
            self.beams += [beam]
            if beam % 2 == 0:
                self.pol2beam['x'] = beam
            else:
                self.pol2beam['y'] = beam
        
        # sort all vcraft files into different polarizations
        all_vfiles=glob.glob(self.vcraft_dr+'/*/*/*.vcraft')
        #nfiles_per_pol = int(len(all_vfiles)/self.npol/self.nant)
        if self.npol == 1:
            self.vfiles = all_vfiles
        elif self.npol == 2:
            vfiles_x = []
            vfiles_y = []
            for vf in all_vfiles:
                if 'beam'+"{:02d}".format(self.pol2beam['x']) in vf:
                    vfiles_x += [vf]
                elif 'beam'+"{:02d}".format(self.pol2beam['y']) in vf:
                    vfiles_y += [vf]
                else:
                    print('ERROR, cannot identify polarization of: '+vf)
            self.vfiles = [vfiles_x,vfiles_y]
            ntotf = len(vfiles_x)+len(vfiles_y)
            if ntotf != len(all_vfiles):
                print('ERROR, file lost in the sorting process. expected: ',len(all_vfiles),', received: ', ntotf)
        else:
            print('ERROR, there should be either 1 or 2 polarizations')

        # get all antenna names for antenna index
        self.an_name = [np.nan]*self.nant
        for n in range(self.nant):
            try:
                self.an_name[n] = self.map_an_name(n)
            except:
                print('Error in maping names for iant',n)
            
        # frequency information for vcraft files
        self.freqs = {}
        freq_offset = -1
        self.nfreq = 0
        for c in range(1,8):
            for f in range(6):
                card_name = 'c'+str(c)+'_f'+str(f)
                with open(self.vcraft_dr+'/'+self.an_name[0]+'/beam'+"{:02d}".format(self.beams[0])+\
                          '/'+self.an_name[0]+'_'+card_name+'.vcraft.hdr', 'r') as fl:
                    for line in fl:
                        if 'FREQS' in line:
                            freqs = (np.array(line.split()[1].split(',')).astype(int) + freq_offset)
                            self.freqs[card_name]=freqs
                            self.nfreq += len(freqs)
        
    def map_an_name(self, an_ind):
        #%% MAP ANTENNA INDEX TO ANTENNA NAME
        # an_name : antenna name in "akNN" format
        with open(self.imfile, 'r') as fl:
            for line in fl:
                if 'TELESCOPE '+str(an_ind)+' NAME: ' in line:
                    an_name = line.split()[-1]
        return an_name
    
    def get_hwdelay(self, an_ind, pol, an_check=0, incards=False):
        #%% DELAY 0: FPGA HARDWARE DELAY
        
        # hwdelay = delays in samples (CLK OFFSET (microsec), frequency dependent)

        search_keyword = "CLOCK COEFF "+str(an_ind)+"/0"

        
        if not incards:
            delays_offset = np.full((self.nfreq),np.nan,dtype=[('frequency',int),('hwdelay', int)])
        else:
            hwdelays = {}

        for c in range(1,8):
            for f in range(6):
                card_name = 'c'+str(c)+'_f'+str(f)
                with open(self.input_dr+'/'+card_name+'/craftfrb.input', 'r') as fl:
                    if an_check: # Quick check to see if antenna index and antenna name match
                        for i in range(self.nant):
                            fl.seek(0)
                            fl.seek(fl.read().find('TELESCOPE NAME '+str(i)),0)
                            input_an_name = str(fl.readline().split()[-1])
                            if input_an_name != self.an_name[i]:
                                print('WARNING: telescope name and index dont match')

                    # read in hw delays

                    vcraft_fn = self.an_name[an_ind]+"_"+card_name+".vcraft"
                    freqs = self.freqs[card_name]
                    # Quick check to compare frequency order with vcraft file
                    for nu in range(len(freqs)):
                        fl.seek(0)
                        fl.seek(fl.read().find('FREQ (MHZ) '+str(nu)),0)
                        data = fl.readline()
                        frequency = int(round(float(data.split()[-1])))
                        #print(frequency)
                        where = np.argwhere(freqs==(frequency-1))[0,0]
                        if where != nu:
                            print("WARNING: the frequency order is different from the vcraft file")


                    fl.seek(0)
                    fl.seek(fl.read().find('TELESCOPE INDEX:    '+str(an_ind)))
                    bool_ypol_yet = False
                    for nu in range(len(freqs)):
                        bool_continue = True
                        if pol == 'x':
                            while bool_continue:
                                line = fl.readline()
                                if 'CLK OFFSET '+str(nu) in line:
                                    data = line
                                    bool_continue = False
                                if "REC BAND 0 POL:     X" in line:
                                    bool_continue = False
                        else:
                            while bool_continue:
                                line = fl.readline()
                                if "REC BAND 0 POL:     X" in line:
                                    bool_ypol_yet = True
                                if bool_ypol_yet:
                                    if 'CLK OFFSET '+str(nu) in line:
                                        data = line
                                        bool_continue = False
                                    if "REC BAND 0 POL:     Y" in line:
                                        bool_continue = False
                        #delays_offset[((c-1)*6+f)*8+nu,0] = freqs[nu]
                        delay_t = float(data.split()[-1])
                        if nu == 0:
                            ref_hwdelay = 8* np.round((delay_t*32/27)/8)
                        hwdelay = 8* int(np.round((delay_t*32/27)/8))
                        if hwdelay != ref_hwdelay:
                            print('WARNING: same fpga has different hwdelay?')
                        if incards:
                            hwdelays[card_name] = hwdelay
                        else:
                            delays_offset[((c-1)*6+f)*8+nu] = (freqs[nu], hwdelay)
        if not incards:
            hwdelays = np.sort(np.copy(delays_offset), order='frequency')
        return hwdelays
    
    
    def get_delay_geo(self,an_ind,SRC=0):
        #%% DELAY 1: GEOMETRIC TIME DELAY
        #SRC = 0 # 0 for pointing centre source, and 1 through to SRC N are the N phase centres

        search_keyword = "SRC "+str(SRC)+" ANT "+str(an_ind)+" DELAY"
        with open(self.imfile, 'r') as fl:
            fl.seek(fl.read().find(search_keyword))
            data=fl.readline()
            delays_geo = float(data.split()[6])/1e6
        return delays_geo
    
    def get_delay_clock(self, an_ind, pol, an_check=0):
        #%% DELAY 2: CLOCK TIME DELAY
        
        # delays_coeff = CLOCK COEFF (microsec) -> unit is seconds


        search_keyword = "CLOCK COEFF "+str(an_ind)+"/0"

        # only for one FPGA
        c=1
        f=0
        card_name = 'c'+str(c)+'_f'+str(f)
        with open(self.input_dr+'/'+card_name+'/craftfrb.input', 'r') as fl:
            if an_check: # Quick check to see if antenna index and antenna name match
                for i in range(self.nant):
                    fl.seek(0)
                    fl.seek(fl.read().find('TELESCOPE NAME '+str(i)),0)
                    input_an_name = str(fl.readline().split()[-1])
                    if input_an_name != self.an_name[i]:
                        print('WARNING: telescope name and index dont match')
            # Read clock delays
            fl.seek(0)
            fl.seek(fl.read().find(search_keyword))
            data=fl.readline()
            delay_clock = float(data.split(':')[-1])

        delay_clock /= 1e6 #convert to seconds
        
        return delay_clock
    
    def get_delay_fring(self, an_ind, pol):
        #%% DELAY 3: FRING FINE TIME DELAY
        # delays_fring : fine time delay measured by FRING
        # sign is different from delays_geo and delays_clock

        if pol == 'x':
            start = "DELAY 1        RATE 1" # Delay 1
            delay_ind = -3 # delay1 location within a line
        else:
            start = "DELAY 2        RATE 2" # Delay 2
            delay_ind = -2 # delay2 location within a line

        with open(self.fringfile, 'r') as fl:

            fl.seek(0)
            fl.seek(fl.read().find(start),0)

            bool_record = False
            while not bool_record:
                data=fl.readline()
                if data.split()[0] == str(an_ind+1):
                    bool_record = True
                    delays_fring = float(data.split()[delay_ind]) # FRING fine time delay

        polfringfile = self.fringfile.replace('delays.sn.txt','xpolfring_xpol.sn')
        if 0: #os.path.exists(polfringfile)==True: # DONT USE THIS NOW. It has been shown that AIPS polarization delay information is not enough. You also need phase offset.
            print('Polarization fring file exists, adding this delay to the fring delay')
            with open(polfringfile, 'r') as fl:
                fl.seek(0)
                fl.seek(fl.read().find(start),0)

                bool_record = False
                while not bool_record:
                    data=fl.readline()
                    if data.split()[0] == str(an_ind+1):
                        bool_record = True
                        delays_polfring = float(data.split()[delay_ind]) # Polarization FRING fine time delay
                        #delays_polfring *= -1
            delays_fring += delays_polfring
            
        return delays_fring
    
    def get_phase_fring(self, an_ind, pol):
        #%% PHASE 1: FRING PHASE
        # phase_fring : phase measured by FRING

        if pol == 'x':
            start = "REAL1          IMAG1" # Phase 1
            delay_ind = 4 # real1 location within a line
        else:
            start = "REAL2          IMAG2" # Phase 2
            delay_ind = 5 # real2 location within a line

        with open(self.fringfile, 'r') as fl:
            fl.seek(0)
            fl.seek(fl.read().find(start),0)

            bool_record = False
            while not bool_record:
                data=fl.readline()
                if data.split()[0] == str(an_ind+1):
                    bool_record = True
                    phase_fring_real = float(data.split()[delay_ind]) # FRING real phase
                    phase_fring_imag = float(data.split()[delay_ind+1]) # FRING imag phase
            phase_fring = phase_fring_real + 1j*phase_fring_imag
            if (abs(phase_fring)-1)>1e-3:
                print("WARNING: amplitude of FRING phase is not 1 but "+str(abs(phase_fring)))
        return phase_fring
    
    def get_phase_selfcal(self, an_ind, pol):
        #%% PHASE 2: SELFCAL PHASE
        # phase_selfcal : phase measured by Selfcal
        # abs_selfcal : amplitude measured by Selfcal

        if pol == 'x':
            start = "REAL1          IMAG1" # Phase 1
            delay_ind = 4 # real1 location within a line
        else:
            start = "REAL2          IMAG2" # Phase 2
            delay_ind = 5 # real2 location within a line

        with open(self.scfile, 'r') as fl:
            fl.seek(0)
            fl.seek(fl.read().find(start),0)
            
            bool_record = False
            while not bool_record:
                data=fl.readline()
                if data.split()[0] == str(an_ind+1):
                    bool_record = True
                    phase_sc_real = float(data.split()[delay_ind]) # selfcal real phase
                    phase_sc_imag = float(data.split()[delay_ind+1]) # selfcal imag phase

            phase_selfcal = phase_sc_real + 1j*phase_sc_imag
        return phase_selfcal
    
    def get_phase_bandpass(self, an_ind, pol):
        #%% PHASE 3: BANDPASS PHASE
        # phase_bandpass : frequency dependent phase
        # abs_bandpass : frequency dependent amplitude

        if pol == 'x':
            start = "REAL 1         IMAG 1" # x pol
            delay_ind = 3 # real1 location within a line
        else:
            start = "REAL 2         IMAG 2" # y pol
            delay_ind = 7 # real2 location within a line
        
        with open(self.bpfile,'r') as fl:
            for line in fl:
                if 'TFDIM11' in line:
                    nfreq = int(line.split()[2])

        with open(self.bpfile, 'r') as fl:
            fl.seek(0)
            fl.seek(fl.read().find(start),0)

            bool_record = False
            phase_bp_real = np.full(nfreq,np.nan)
            phase_bp_imag = np.full(nfreq,np.nan)
            while not bool_record:
                data=fl.readline()
                if data.split()[0] == str(an_ind+1):
                    bool_record = True
                    for chan in range(nfreq):
                        phase_bp_real[chan] = float(data.split()[delay_ind]) # bandpass real phase
                        phase_bp_imag[chan] = float(data.split()[delay_ind+1]) # bandpass imag phase
                        data = fl.readline()

            phase_bandpass = phase_bp_real + 1j*phase_bp_imag
        return phase_bandpass



