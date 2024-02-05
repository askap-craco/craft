#!/usr/bin/env python
"""
Vismetaata data holder class

Copyright (C) CSIRO 2022
"""
import sys
import logging
from abc import ABC
from craft.craco import ant2bl,get_max_uv
from craft.freq_config import FrequencyConfig
from astropy.time import Time
from astropy.coordinates import SkyCoord

log = logging.getLogger(__name__)

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

class VisMetadata:
    '''
    Data class for vis metadata

    You can give one of these to CRACO plan and it will make new plan for you

    :baselines: dictionary of baselines. Key=float, baseline ID (UVFITS) value=Something with {UU','VV','WW' in seconds)
    :freq_config: frequency configuration - see freq config object
    :target_name:str
    :target_skycoord: target sky coordinate
    :beamid: Beam id (int)
    :tstart: Centre of tstart start time of first integration
    :tsamp: Quantity with units of time
    :tuv: Time at which UV data was valid - different from tstart possibly
    '''
    def __init__(self,
                 baselines:dict,
                 freq_config:FrequencyConfig,
                 target_name:str,
                 target_skycoord:SkyCoord,
                 beamid:int,
                 tstart:Time,
                 tsamp,
                 tuv:Time):
        self._baselines = baselines
        self._freq_config = freq_config
        self._target_name = target_name
        self._target_skycoord = target_skycoord
        self._beamid = beamid
        self._tstart = tstart
        self._tsamp = tsamp
        self._flag_ants = []
        self._tuv = tuv

    def set_flagants(self, flag_ants):
        '''
        Set flag antennas as a list of 1-based integers
        '''
        self._flag_ants = flag_ants

    def get_max_uv(self) -> float:
        ''' Return umax, vmax'''
        fmax = self.channel_frequencies.max()
        baselines = self.baselines
        maxuv = get_max_uv(baselines, fmax)
        return maxuv


    @property
    def baselines(self) -> dict:
        '''
        Returns dictionary by blid of baselines, containing UU,VV,WW in seconds
        Like a uv fits file
        '''
        return self._baselines

    @property
    def nbl(self) -> int:
        return len(self.baselines)

    @property
    def channel_frequencies(self):
        '''
        Return np array of channel frequencies in Hz.
        '''
        return self.freq_config.channel_frequencies*1e6

    @property
    def freq_config(self) -> FrequencyConfig:
        '''
        Returns the frequency configuration object
        '''
        return self._freq_config

    @property
    def target_name(self) -> str:
        '''
        String target name
        '''
        return self._target_name

    @property
    def target_skycoord(self) -> SkyCoord:
        '''
        Return target as SkyCoord
        '''
        return self._target_skycoord

    @property
    def beamid(self) -> int:
        return self._beamid

    @property
    def tstart(self) -> Time:
        return self._tstart

    @property
    def tsamp(self):
        '''
        Return sample interval in .... seconds?
        '''
        return self._tsamp

    @property
    def tuv(self) -> Time:
        return self._tuv

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    

if __name__ == '__main__':
    _main()
