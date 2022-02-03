#!/usr/bin/env python
"""
Classes for extracting information from CRAFT real-time pipeilne.

Copyright (C) CSIRO 2019
"""
import numpy as np
import os
import sys
import logging
from . import crafthdr
import re
import glob
from . import sigproc
from . import dada

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

fred_file_re = re.compile('(.*).dada')

def type_of_file(p):
    f = os.path.basename(p)
    m = fred_file_re.match(f)
    if m is None:
        raise ValueError('Path {} is not a valid FREDDA output'.format(p))

    return g[0]

class FreddaRescaleBlock(dict):
    def __init__(self, rsdata, blkid):

        for df in rsdata.dada_files:
            name = df.data_name
            data = df[blkid]
            # most data types are size [1,nant,nbeams,nf] except dm0 and dm0count which is [1,nant,nbeams,nt] and dm0stats which is [1,nant,nbeams,4]
            # in all cases, we can remove the first dimension so as not to
            # confuse the user
            
            assert data.shape[0] == 1
            data = data[0,rsdata.antenna_idxs,:,:]
            self[name] = data

        self.mjd = rsdata.tstart + blkid*rsdata.tsamp/86400.
        self.rsdata = rsdata
        self.blkid = blkid

class FreddaRescaleData(object):
    def __init__(self, path, exclude_ants=None):
        self.path = path

        dada_file_paths = glob.glob(os.path.join(self.path, '*.dada'))
        self.dada_files = [dada.DadaFile(f) for f in dada_file_paths]
        hdr = self.dada_files[0].hdr
        self.hdr = hdr
        
        f1 = float(hdr.get_value('SOURCE_FCH1'))
        foff = float(hdr.get_value('SOURCE_FOFF'))
        nant = int(hdr.get_value('SOURCE_NANTS'))
        nchan = int(hdr.get_value('NF'))
        self.freqs = np.arange(nchan)*foff + f1
        self.nbeams_per_antenna = int(self.hdr.get_value('NBEAMS_PER_ANTENNA'))
        self.nbeams = int(hdr.get_value('NBEAMS_OUT'))
        self.nbeams_per_antenna = int(self.hdr.get_value('NBEAMS_PER_ANTENNA'))
        self.npol = int(self.hdr.get_value('NPOLS_IN'))
        self.nbteams_in_total = int(self.hdr.get_value('NBEAMS_IN_TOTAL'))
        self.nt = int(self.hdr.get_value('NT'))
        self.nd = int(self.hdr.get_value('ND'))
        self.tstart = float(self.hdr.get_value('SOURCE_TSTART'))
        self.tsamp = float(self.hdr.get_value('TSAMP')) # tsamp fo rrescaling different from source tsamp
        self.antennas = self.hdr.get_value('SOURCE_ANTENNA_NAME','')
        if self.antennas.strip() == '':
            antennas = ['ia{:02d}'.format(ia) for ia in range(nant)]
        else:
            antennas = self.antennas.split(',')

        if exclude_ants is None:
            self.antennas = antennas
            self.antenna_idxs = slice(None) # select everything
        else:
            self.antennas = []
            self.antenna_idxs = []
            for ia, a in enumerate(antennas):
                antid = int(a[2:])
                if antid not in exclude_ants:
                    self.antennas.append(a)
                    self.antenna_idxs.append(ia)
                
    @property
    def nblocks(self):
        return min(list(map(len, self.dada_files)))

    def get_block(self, blkid):
        return FreddaRescaleBlock(self, blkid)

    def __getitem__(self, blkid):
        return self.get_block(blkid)

    def blocks(self, step=1):
        assert step >= 1
        blkid = 0
        while blkid < self.nblocks:
            yield self[blkid]
            blkid += step


class DataDir(object):
    def __init__(self, rootobj, path):
        self.rootobj = rootobj
        self.path = path
        assert os.path.isdir(self.path)
        self.name = os.path.split(self.path)[-1].replace('/','')

    def _pathjoin(self, *args):
        return os.path.join(self.path, *args)
    
    def _pathglob(self, pattern, isdir=False):
        lst = glob.glob(self._pathjoin(pattern))
        if isdir:
            lst = list(filter(os.path.isdir, lst))

        return lst

    def __str__(self):
        return 'DataDir:'+str(self.path)

    __repr__ = __str__

class StartDir(DataDir):

    def __init__(self, *args, **kwargs):
        DataDir.__init__(self, *args, **kwargs)
        self._hdr = None # lazy
        
    @property
    def antname(self):
        return self.rootobj.name

    @property
    def scanname(self):
        return self.rootobj.rootobj.name

    @property
    def sbname(self):
        return self.rootobj.rootobj.rootobj.name

    @property
    def start_fullname(self):
        return os.path.join(self.sbname, self.scanname, self.antname, self.name)

    @property
    def header(self):
        if self._hdr is not None:
            return self._hdr
        
        hdrname = self._pathglob('{}*.hdr'.format(self.antname))
        if len(hdrname) == 0:
            #raise ValueError('No valid headers in {}'.format(self.path))
            hdr = None
        elif len(hdrname) > 1:
            raise ValueError('Too many headers in {} = {}'.format(self.path, ' '.join(hdrname)))
        else:
            hdr = crafthdr.DadaHeader.fromfile(hdrname[0])

        self._hdr = hdr
            
        return hdr

    @property
    def filterbank_paths(self):
        filnames = self._pathglob('*.[0-9][0-9].fil')
        return filnames

    @property
    def num_beams(self):
        return len(self.filterbank_names)

    @property
    def filterbanks(self):
        f = [sigproc.SigprocFile(f) for f in self.filterbank_paths]
        return f

    def open_filterbank(self, beamid):
        '''
        Opens a filterbank for the given beamid
        '''
        fname = self._pathglob('*.{:02d}.fil'.format(beamid))
        if len(fname) == 0:
            raise ValueError('No beamid = {} in path {}'.format(beamid, self.path))

        assert len(fname) == 1, 'Multiple filterbanks in path {}'.format(self.path)

        s = sigproc.SigprocFile(fname[0])
        return s

    @property
    def fredda_candfile_name(self):
        '''
        Returns a path to the fredda candidate file name if it exists. Otherwise None
        '''
        
        f = self._pathjoin('fredda.cand')
        if os.path.exists(f):
            return f
        else:
            return None

    def _load_fredda_output(sel, name, epoch):
        fname = self._pathjoin('{}_e{:d}.dat'.format(name, epoch))
        d = load4d(fname)

        return d

    def list_fredda_output_paths(self):
        '''
        Returns a list of all the paths to fredda outputs.
        If epoch=None, then it will return all epochs
        Otherwise, specific epoch as an int
        It there are no fredda outputs, it will return an empty list
        '''
        paths = self._pathglob('*.dada')

        return paths

    def list_fredda_output_types(self, epoch=None):
        '''
        Returns a set of the types of FREDDA output available a the given epoch.
        e.g. ['mean','std','kurt','inbuf','fdmt', 'dm0','inbuf','nsamps']
        if epoch is None, it returns tpes for all epochs.
        '''
        
        paths = self.list_fredda_output_paths(epoch)
        fred_types = set(type_of_file(p) for p in paths)

        return fred_types

    def fredda_data(self):
        return FreddaRescaleData(self.path)

class ScanDir(DataDir):
    
    @property
    def antenna_paths(self):
        return self._pathglob('[ak|co]*', isdir=True)
    
    @property
    def antenna_names(self):
        if self._ants is None:
            self._ants = [os.paths.split(d)[-1] for d in self.antenna_paths]
        return self._ants

    @property
    def antenna_dirs(self):
        adirs = [AntennaDir(self, d) for d in self.antenna_paths]
        return adirs

    @property
    def start_names(self):
        '''
        Returns a set of all the start names. Not all antennas will have all starts
        '''
        all_starts = []
        for a in self.antenna_dirs:
            all_starts.extend(a.start_names)

        return set(all_starts)

    def get_antennas_with_start(self, start_name):
        ''' 
        Returns all the antennas with the given start
        '''
        ants = [d for d in self.antenna_dirs if start_name in d.start_names]
        return ants
                

    @property
    def is_ics(self):
        '''
        Returns True if this start has an incoherent sum ('ICS') directory
        '''
        has_ics =  os.path.isdir(self._pathjoin('ICS'))
        return has_ics

    def __str__(self):
        return 'ScanDir: {} ICS={}'.format(self.path, self.is_ics)


class AntennaDir(DataDir):

    @property
    def start_paths(self):
        paths = self._pathglob('C*', isdir=True)
        return paths

    @property
    def start_dirs(self):
        sta = [StartDir(self, d) for d in self.start_paths]
        return sta

    @property
    def start_names(self):
        return list(map(os.path.basename, self.start_paths))

    def get_start(self, start_name):
        return StartDir(self, self._pathjoin(start_name))

    @property
    def is_open(self):
        '''
        Returns true if the scan dir is currently open
        '''
        return os.path.exists(self._pathjoin('SCAN_OPEN'))

    @property
    def is_closed(self):
        '''
        Returns True if the scan dir is closed.
        Probably more secure to run not self.is_open rather than self.is_closed
        '''
        return os.path.exists(self._pathjoin('SCAN_CLOSED'))

    def __str__(self):
        s = 'ScanDir: '
        s += self.path
        if self.is_open:
            s += 'OPEN'
            
        if self.is_closed:
            s += 'CLOSED'

        return s

class SchedblockDir(DataDir):

    @property
    def scan_paths(self):
        '''
        Returns a list of fully qualified scan paths as strings
        e.g. /some/data/SB123/2019010203040405')
        '''
        scannames = self._pathglob('2*', isdir=True)
        return scannames

    @property
    def scan_names(self):
        '''
        Retuns a list of scan names as strings. Only the scan npart of the name, not the other bits
        e.g. ['201901020304050607')
        '''
        scannames = list(map(os.path.basename, self.scan_paths))
        return scannames
                    

    @property
    def scan_dirs(self):
        thescans =  [ScanDir(self, d) for d in self.scan_paths]
        return thescans

class SearchDir(DataDir):
    @property
    def sbpaths(self):
        sbnames = self._pathglob('SB*', isdir=True)
        return sbnames

    @property
    def sbnames(self):
        return list(map(os.path.basename, self.sbpaths))

    @property
    def sbdirs(self):
        sblist = [SchedblockDir(self, d) for d in self.sbpaths]
        return sblist

    def get_schedblock(self, sbname):
        return SchedblockDir(self, self._pathjoin(sbname))
    

class RtData(object):
    def __init__(self, data_dir=None, voltage_dir=None):
        self.data_dir = data_dir
        if self.data_dir is None:
            self.data_dir = os.environ['CRAFT_DATA_DIR']

        if not os.path.isdir(self.data_dir):
            raise ValueError('Data dir is not a directory {}'.format(self.data_dir))

        self.searchdata = SearchDir(self, self.data_dir)
        
        self.voltage_dir = voltage_dir
        if self.voltage_dir is None:
            self.voltage_dir = os.environ.get('CRAFT_VOLTAGE_DIR', None)

        self.voltages = None
        if self.voltage_dir is not None:
            self.voltages = DataDir(self, self.voltage_dir)

    def __str__(self):
        return 'RtData: search={} voltage={}'.format(self.searchdata, self.voltages)


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('commands', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    command = values.commands[0]

    d = RtData()
    if command == 'listsb':
        for sb in d.searchdata.sbnames:
            print(sb)
    if command == 'showsb':
        sb = d.searchdata.get_schedblock(values.commands[1])
        print(sb)
        for scan in sb.scan_dirs:
            print((' '*1+ str(scan)))
            for ant in scan.antenna_dirs:
                print((' '*2 + str(ant)))
                for start in ant.start_dirs:
                    print((' '*3 + str(start)))
                    for fil in start.filterbanks:
                        print((' '*4 + str(fil)))
            for start in scan.start_names:
                print(start)
                ants = scan.get_antennas_with_start(start)
                antstarts = [a.get_start(start) for a in ants]
                for antstart in antstarts:
                    print(antstart.start_fullname)
                    print(antstart.header)
                    print(antstart.fredda_candfile_name)
                    for fredepoch in antstart.list_fredda_epochs:
                        print(' '*5, str(fredepoch))
                    

                    for f in antstart.filterbank_paths:
                        print(' '*5, str(f))


                        
    
            
        
    

if __name__ == '__main__':
    _main()
