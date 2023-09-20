#!/usr/bin/env python
"""
Plans a CRACO scan

Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import pickle 
import warnings
from astropy.wcs import WCS
from astropy import units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, Angle

from craft.craco_wcs import CracoWCS

from . import uvfits
from . import craco_kernels
from . import craco
from . import fdmt
from .fdmt_plan import *
from .cmdline import strrange
from collections import namedtuple

log = logging.getLogger(__name__)

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

BMAX = 6e3; # max baseline meters
FREQMAX = 1.8e9 # max frequency Hz
C = 3e8 # m/s
LAMMIN = C/FREQMAX # smallest wavelength
UVMAX = BMAX/LAMMIN # largest UV point

def splitn(s, delim, n=2):
    '''
    Splits by the  delimiter. If there's only one value, repeat it by n
    '''
    bits = s.split(delim)
    if len(bits) == 1:
        bits = bits*n

    return bits

class ImageParams1d:
    def __init__(self, bmax:float, npix:int, os:float, fov:Angle, imgpix=None):
        '''
        Calculates image parameters for a given axis
        normally: Tries to achieve the desired Field of View. If it can't it produces an 
        FoV that has the desired oversampling. 
        If image pixel cell size (in arcseconds) is specified, then it just uses that
    
        :bmax: max baseline length in units of lambda
        :npix: image size pixels
        :os: desired minimum oversampling of pixels per synthesized beam
        :fov: desired FoV as an Angle
        :imgpix: Image pixel size (Angle). Default=None
        '''

        synth_beam = Angle(1/bmax, u.radian) # Synthesized beam size in units of lambda
    
        if imgpix is not None: # use image pixel size directly
            img_pix_size = imgpix
            actual_os = synth_beam.deg / img_pix_size.deg
            if actual_os < os:
                raise ValueError('Given imgpix cant achieve desired oversampling')
        else:
            actual_os = synth_beam.deg * npix / fov.deg
            if actual_os < os: # could not achieve desired FoV
                actual_os = os # set the actual oversampling
                
            img_pix_size = synth_beam / actual_os

        actual_fov = img_pix_size * npix
        uv_cell_size = 1/actual_fov

        self.fov = actual_fov
        self.requested_fov = fov
        
        self.os = actual_os
        self.requested_os = os

        self.requested_img_pix_size = imgpix
        self.img_pix_size = img_pix_size
        self.uv_cell_size = uv_cell_size
        self.synth_beam = synth_beam
        self.bmax = bmax

    def __str__(self):
        s = f'ImageParms1D fov={self.fov.deg:0.2f}({self.requested_fov.deg:0.2f})deg os={self.os:0.2f}({self.requested_os:0.2f}) imgpix={self.img_pix_size.arcsec:0.1f} arcsec synth beam={self.synth_beam.arcsec:0.1f} arcsec uvcell={self.uv_cell_size.value:0.1f} lambda'
        return s

class ImageParams2d:
    def __init__(self, uvmax, npix, os_str, fov_str):
        self.uvmax = uvmax
        self.npix = npix
        requested_os = list(map(float, splitn(os_str, ',', 2)))
        requested_fov = list(map(Angle, splitn(fov_str, ',', 2)))
        self.lparams = ImageParams1d(uvmax[0], npix, requested_os[0], requested_fov[0])
        self.mparams = ImageParams1d(uvmax[1], npix, requested_os[1], requested_fov[1])

    @property
    def lmcell(self):
        '''
        Returns 2-Angle with units in arcsec
        '''
        return Angle([self.lparams.img_pix_size, self.mparams.img_pix_size], u.arcsec)

    @property
    def uvcell(self):
        '''
        Returns 2-tuple as raw values in lambda
        '''
        return (self.lparams.uv_cell_size.value, self.mparams.uv_cell_size.value)

    def __str__(self):
        s = f'ImageParams2d lparams={self.lparams} mparams={self.mparams}'
        return s
    
    
def dump_plan(plan, pickle_fname):
    filehandler = open(pickle_fname, 'wb') 
    pickle.dump(plan, filehandler)
    filehandler.close()

def load_plan(pickle_fname):        
    filehandler = open(pickle_fname, 'rb')
    plan = pickle.load(filehandler)
    filehandler.close()
    
    return plan
            
def make_ddreader_configs(plan):
    '''
    Make the ddgrid reader config
    See maek_sample_configs in ddgrid_catchtest.cpp
    '''

    configs = np.zeros((plan.nuvrest_max, 2), dtype=np.uint16)
    for irun, run in enumerate(plan.fdmt_plan.runs):
        # now convert to integers
        if run is None: # I'm not sure whether 0 makes sense but we'll work it out for now
            configs[irun, 0] = 0
            configs[irun, 1] = 0
        else:
            configs[irun, 0] = fdmt.cff_to_word(run.idm_cff)
            configs[irun, 1] = fdmt.cff_to_word(run.offset_cff)

    return configs


def fftshift_coordinate(c, npix):
    '''
    Returns the FFT shifted version of the given coordinate with respect to an FFT that has size npix

    >>> fftshift_coordinate(0, 256)
    128
    >>> fftshift_coordinate(128, 256)
    0
    >>> fftshift_coordinate(1, 256)
    129
    >>> fftshift_coordinate(129, 256)
    1
    '''
    assert c >= 0
    pivot = npix // 2
    assert npix % 2 == 0
    if c < pivot:
        cout = c + pivot
    else:
        cout = c - pivot

    assert 0 <= cout < npix
    
    return cout

fftshift_coordinates = np.vectorize(fftshift_coordinate)

def get_uvcells(baselines, uvcell, freqs, Npix, plot=False, fftshift=True, transpose=False):
    '''
    fftshift=True means apply similar transform to numpy fft shift so that the DC bin is the first elemetn
    '''
    uvcells = []

    ucell, vcell = uvcell

    if plot:
        grid = np.zeros((Npix, Npix))

    for blid, bldata in list(baselines.items()):
        #UU, VV WW are in seconds
        ulam = bldata['UU'] * freqs
        vlam = bldata['VV'] * freqs
        a1,a2 = craco.bl2ant(blid)

        if np.any((np.abs(ulam) >= UVMAX) | (np.abs(vlam >=UVMAX))):
            warnings.warn(f'Maximum UV is larger than ASKAP for {a1}-{a2}. UVMAX={UVMAX} ulam={ulam.max()} vlam={vlam.max()}')

        pix_offset = float(Npix)/2.0

        # OK - so the UV plane is looking down on the telescope
        # with east to the right, but we wan tot make the image
        # with east to the left as though we're looking up at it.
        # so we apply a minus sign to U to make it negative.
        # hopfully this works
        upix = np.round(-ulam/ucell + pix_offset).astype(int)
        vpix = np.round(vlam/vcell + pix_offset).astype(int)
        if np.any((upix < 0) | (upix >= Npix) | (vpix < 0) | (vpix >= Npix)):
            warnings.warn('Pixel coordinates out of range')
            raise ValueError(f'Pixel coordinates out of range, upix = {upix}, vpix = {vpix}')

        if fftshift:
            upix = fftshift_coordinates(upix, Npix)
            vpix = fftshift_coordinates(vpix, Npix)

        if transpose:
            (vpix, upix) = (upix, vpix)

        if plot:
            #pylab.plot(ulam/1e3, vlam/1e3)
            pylab.plot(ulam/ucell + pix_offset, vlam/vcell+pix_offset)

        uvpos = list(zip(upix, vpix))
        for istart, iend in craco.runidxs(uvpos):
            uvpix = uvpos[istart]
            assert uvpos[istart] == uvpos[iend]
            assert uvpix[0] < Npix
            assert uvpix[1] < Npix
            b = craco.BaselineCell(blid, uvpix, istart, iend, freqs[istart:iend+1], Npix)
            if uvpix[0] == 0 or uvpix[1] == 0:
                warnings.warn(f'Cannot grid things on U=0 or V=0 blid={blid} {a1}-{a2} uvpix={uvpix}')
            else:
                uvcells.append(b)
                #print(b.uvpix_lower, b.uvpix_upper, ulam[0], vlam[0], ulam[0]/ucell, vlam[0]/vcell)#, ulam*freqs[0]/ucell, vlam*freqs[0]/vcell)
                #print(b)

                
            if plot:
                grid[uvpix[1],uvpix[0]] += len(b.freqs)


    if plot:
        pylab.imshow(grid)
        pylab.show()

    uvcells = sorted(uvcells, key=lambda b:b.upper_idx)
    return uvcells

class AddInstruction(object):
    def __init__(self, plan, target_slot, cell_coords, uvpix):
        self.plan = plan
        self.target_slot = target_slot
        self.cell_coords = cell_coords
        self.uvpix = uvpix
        self.shift = False

    @property
    def shift_flag(self):
        return 1 if self.shift else 0

    @property
    def uvidx(self):
        irun, icell = self.cell_coords
        c = icell + self.plan.nuvwide*irun
        return c

    def __str__(self):
        irun, icell = self.cell_coords
        cell = self.plan.fdmt_plan.get_cell(self.cell_coords)
        return 'add {self.cell_coords} which is {cell} to slot {self.target_slot} and shift={self.shift}'.format(self=self, cell=cell)

    __repr__ = __str__

def calc_grid_luts(plan, upper=True):
    uvcells = plan.uvcells

    # Calculate upper or lower pixel positions but always sort in the same raster order left-to-right, top-to-bottom
    if upper:
        uvpix_list = [cell.uvpix_upper for cell in uvcells]
        sorter = lambda c: craco.triangular_index(c[0],c[1], plan.npix)
    else:
        uvpix_list = [cell.uvpix_lower for cell in uvcells if cell.uvpix[0] != cell.uvpix[1]] # don't include diagonal
        sorter = lambda c: craco.triangular_index(c[1],c[0], plan.npix, raster='yx') # triangular_index requires upper_triangular coordinates
        
    unique_uvcoords = set(uvpix_list)
    # sorts in raster order

    unique_uvcoords = sorted(unique_uvcoords, key=sorter)
    ncoord = len(unique_uvcoords)
    log.info('Got %d unique UV coords. Upper=%s', ncoord, upper)

    # check
    check = False
    if check:
        grid = np.zeros((plan.npix, plan.npix))
        for i, (u, v) in enumerate(unique_uvcoords):
            grid[v,u] = i+1

        pylab.imshow(np.ma.masked_equal(grid, 0), interpolation='none', origin='lower')
        pylab.title(f'Checking calc_grid_luts upper={upper}')
        pylab.show()
    
    fplan = plan.fdmt_plan
    fruns = fplan.runs
    remaining_fdmt_cells = [] # this is an array we keep for self-testing to make sure we used everything
    for run in fruns:
        if run is not None:
            remaining_fdmt_cells.extend(run.defined_cells)

    ngridreg = plan.ngridreg
    all_instructions = []

    nwrites = (ncoord + ngridreg - 1) // ngridreg # number of writes
    assert nwrites <= 4096*2, 'Not enough clocks to write data in! nwrites={}'.format(nwrites)
    log.info('Need to write %d groups of %d register to pad function', nwrites, ngridreg)

    for iwrite in range(nwrites):
        # Add available cells
        n = min(ngridreg, ncoord - iwrite*ngridreg)
        

        slot_uv_coords = unique_uvcoords[iwrite*ngridreg:iwrite*ngridreg + n]
        log.debug('n=%s ncoord=%s', n, len(slot_uv_coords))
        instructions = []

        for islot, slot_uv_coord in enumerate(slot_uv_coords):
            # The FDMT has all its data in the upper hermetian - so we always use upper hermetian coordinates
            upper_slot_uv_coord = craco.make_upper(slot_uv_coord, plan.npix)
            log.debug('Considering islot=%d slot=%s upper=%s', islot, slot_uv_coord, upper_slot_uv_coord)

            # Find irun and icell 
            for cell_coord in fplan.find_uv(upper_slot_uv_coord):
                # cell_coord is the location in the FDMT output buffer
                cell = fplan.get_cell(cell_coord) # this is the actual baseline cell. Don't use it other than to remove it form our check buffer
                log.debug('removing cell %s %s', cell, cell_coord)
                remaining_fdmt_cells.remove(cell)

                inst = AddInstruction(plan, islot, cell_coord, slot_uv_coord)
                instructions.append(inst)

        # Grid reading reads 2 UV/clk, so we need to add a dummy add instruction
        if len(instructions) % 2 == 1: # if odd
            instructions.append(AddInstruction(plan, ngridreg-1, fplan.zero_cell, (-1,-1)))

        # See shift marker for last 2 instruction, as we read 2UV/clk and I don't know which mattes
        #instructions[-2].shift = True
        assert len(instructions) > 0
        instructions[-1].shift = True

        all_instructions.extend(instructions)
        log.debug('Added %d instructions with %d available', len(instructions), len(remaining_fdmt_cells))

    # Check instructions

    assert all_instructions[-1].shift == True
    #assert all_instructions[-2].shift == True
    if upper:
        assert len(remaining_fdmt_cells) == 0, f'Leftover cells {remaining_fdmt_cells}'

    num_shifts =  sum([inst.shift == True for inst in all_instructions])

    unique_uvidxs = set([inst.uvidx for inst in all_instructions])
    #assert len(unique_uvidxs) == len(uvcells), 'Got {} unique UV indexces != len(uvcels) = {}'.format(len(unique_uvidxs), len(uvcells))

    assert num_shifts == nwrites, 'num_shifts={num_shifts} != nwrites={nwrites}'.format(**locals())

    all_uvpixs = set([inst.uvpix for inst in all_instructions])
    #assert len(all_uvpixs) == len(all_instructions), 'Got %d unique uvpixels and %d instructions' % (len(all_uvpixs), len(all_instructions))

    return all_instructions

def get_pad_input_registers(instr, ssr=16):
    ''' 
    Groups the inputs by the .shift attribute adn returns the value of of the uvpix attribute
    And returns the groups - only the uvpix part of the register
    '''

    #curr = [None for i in xrange(ssr)]
    curr = [None for i in range(ssr)]
    for i, instruction in enumerate(instr):
        if instruction.uvpix == (-1,-1):
            # add zero instruction. just pass
            pass
        else:
            assert curr[instruction.target_slot] is None or curr[instruction.target_slot] == instruction.uvpix, 'Invalid instruction i={} curr={} instruction={}'.format(i, curr, instruction)
            curr[instruction.target_slot] = instruction.uvpix
            
        if instruction.shift:
            yield curr
            #curr = [None for i in xrange(ssr)]
            curr = [None for i in range(ssr)]

def calc_pad_lut(plan, ssr=16):
    '''
    Calculates the padding lookup table

    @param plan PipelinePlan
    @param ssr super-sample rate for FFT. I.e. how many pixels per clock we need to produce after padding
    @returns (upper_idxs, upper_shifts, lower_idxs, lower_shifts) tuple
    '''
    upper_inst = plan.upper_instructions
    lower_inst  = plan.lower_instructions
    upper_inputs = get_pad_input_registers(upper_inst)
    lower_inputs = get_pad_input_registers(lower_inst)
    uvpix_set = set([cell.uvpix_upper for cell in plan.uvcells])

    lut = []
    upper_registers = []
    lower_registers = []

    # add 2 sets of inputs to upper and lower registers

    upper_registers.extend(next(upper_inputs))
    upper_registers.extend(next(upper_inputs))
                                            
    lower_registers.extend(next(lower_inputs))
    lower_registers.extend(next(lower_inputs))

    log.debug('Upper registers: %s', upper_registers)
    log.debug('Lower registers: %s', lower_registers)
    
    upper_shifts = []
    lower_shifts = []
    upper_idxs = []
    lower_idxs = []

    npix = plan.npix

    # raster scan in u direction fastest, then v
    # top left pixel is 0,0
    for v in range(npix):
        for ublk in range(npix//ssr):
            lower_shift = False
            upper_shift = False
            for iu in range(ssr):
                u = iu + ublk*ssr
                uv = (u,v)
                if u >= v: # upper hermetian
                    if (u,v) in uvpix_set:
                        # it better be in the current register
                        i = upper_registers.index(uv)
                        # For debugging
                        upper_registers[i] = 'USED'
                        upper_shift |= i >= ssr - 1
                    else:
                        i = -1 # which means pad with 0

                    upper_idxs.append((u,v,i))
                else: # lower hermeitan
                    if craco.make_upper(uv, npix) in uvpix_set:
                        i = lower_registers.index(uv)
                        lower_shift |= i >= ssr - 1
                    else:
                        i = -1
                        
                    lower_idxs.append((u,v,i))


            if upper_shift:
                upper_registers[:ssr] = upper_registers[ssr:]
                try:
                    upper_registers[ssr:] = next(upper_inputs)
                except StopIteration:
                    upper_finished = True


            if lower_shift:
                lower_registers[:ssr] = lower_registers[ssr:]
                try:
                    lower_registers[ssr:] = next(lower_inputs)
                except StopIteration:
                    lower_finished = True

            upper_shifts.append(upper_shift)
            lower_shifts.append(lower_shift)


    assert upper_finished, 'Upper inputs not empty'
    assert lower_finished, 'Lower inputs not empty'

    assert len(upper_shifts) == npix*npix/ssr

    return (upper_idxs, upper_shifts, lower_idxs, lower_shifts)


class PipelinePlan(object):
    def __init__(self, f, values=None, dms=None, prev_plan=None):
        '''
        Creates a pipeline plan
        :f: object that contains lots of juicy info like baselines, frequency, tsamp etc. See uvfits.py and search_pipeline_sink:Adapter
        :values: command line arguments
        :dms: list of dispersion measures
        :prev_plan: previous plan in time - todo check this makes sense, but mainly use to work out new FdmtPlan based on old FDMT plan
        
        '''

        if values is None:
            print('Creating default values')
            self.values = get_parser().parse_args()
        elif isinstance(values, str):
            print(f'parsing values {values}')
            self.values = get_parser().parse_args(values.split())
        else:
            self.values = values

        self.prev_plan = prev_plan
        
        values = self.values
        if values.flag_ants:
            f.set_flagants(values.flag_ants)
            
        try:
            beamid = f.beamid
        except:
            log.info('Unknown beamid')
            beamid = -1

        self.beamid = beamid
        
        log.info('making Plan values=%s prev plan:%s', self.values, self.prev_plan)
        self.__tsamp = f.tsamp
        umax, vmax = f.get_max_uv()
        uvmax = (umax,vmax)
        baselines = f.baselines
        self.baselines = baselines
        # List of basleine IDs sorted
        # THis is the ordering of baselines once you've done bl2array
        self.baseline_order = sorted(self.baselines.keys())

        nbl = len(baselines)
        freqs = f.channel_frequencies
        self.target_name = f.target_name

        # Cant handle inverted bands - this is assumed all over the code. It's boring
        assert freqs.min() == freqs[0]
        assert freqs.max() == freqs[-1]
        Npix = self.values.npix
        self.npix = Npix
        self.image_params = ImageParams2d(uvmax, self.npix, self.values.os, self.values.fov)
        self.lmcell = self.image_params.lmcell
        self.uvcell = self.image_params.uvcell

        fmax = freqs.max()
        foff = freqs[1] - freqs[0]
        lambdamin = 3e8/fmax
        umax_km = umax*lambdamin/1e3
        vmax_km = vmax*lambdamin/1e3

        # Could get RA/DeC from fits table, or header. Header is easier, but maybe less correct

        target_skycoord = f.target_skycoord
        self.phase_center = target_skycoord
        self.ra = target_skycoord.ra
        self.dec = target_skycoord.dec
        self.tstart = f.tstart
        craco_wcs = CracoWCS.from_plan(self)
        self.craco_wcs = craco_wcs
        self.wcs = craco_wcs.wcs2 # The 2D WCS for images

        log.info('Nbl=%d Fch1=%f foff=%f nchan=%d lambdamin=%f uvmax=%slambda max baseline=%skm wcs=%s image_params=%s',
                 nbl, freqs[0], foff, len(freqs), lambdamin, (umax, vmax), (umax_km, vmax_km), self.wcs,
                 self.image_params)
        
        uvcells = get_uvcells(baselines, self.uvcell, freqs, Npix, values.show)
        if umax > UVMAX or vmax > UVMAX:
            raise ValueError(f'Maximum basline larger than ASKAP umax={umax} vmax={vmax} UVMAX={UVMAX}')
        if umax < 00/0.2 or vmax < 100/0.2:
            raise ValueError(f'Maximum baseline is unreasonably small. wrong units? umax={umax} vmax={vmax} in km={umax_km}/{vmax_km}')
                  
        log.info('Got Ncells=%d uvcells', len(uvcells))
        d = np.array([(v.a1, v.a2, v.uvpix[0], v.uvpix[1], v.chan_start, v.chan_end) for v in uvcells], dtype=np.int32)
        if self.values.uv is not None:
            np.savetxt(self.values.uv+'.uvgrid.txt', d, fmt='%d',  header='ant1, ant2, u(pix), v(pix), chan1, chan2')

        self.uvcells = uvcells
        self.nt = self.values.nt
        self.ncu = 4 # hard coded
        self.nchunk_time = self.values.nt // (2*self.ncu)
        self.freqs = freqs
        self.nbox = self.values.nbox
        self.boxcar_weight = self.values.boxcar_weight
        self.nuvwide = self.values.nuvwide
        self.nuvmax  = self.values.nuvmax
        assert self.nuvmax % self.nuvwide == 0
        self.nuvrest_max = self.nuvmax // self.nuvwide
        self.ncin  = self.values.ncin
        self.ndout = self.values.ndout
        self.foff = foff
        self.dtype = np.complex64 # working data type
        self.nbl = nbl
        self.baseline_shape = (self.nbl, self.nf, self.nt)

        self.fdmt_scale = self.values.fdmt_scale
        self.fft_scale  = self.values.fft_scale
        
        self.fft_ssr = 16 # number of FFT pixels per clock - "super sample rate"
        self.ngridreg = 16 # number of grid registers to do

        # calculate DMs and DD grid reader LUT
        # DMS is in units of samples
        # you could load from a file or a more difficult specification from the CMDLINE
        # For now we just do every DM up to nd
        assert self.values.ndm <= self.values.max_ndm, 'Requested ND is > larger than MAX_NDM'
        if dms is None:
            dms = np.arange(self.values.ndm, dtype=np.uint32)

        self.__set_dms(dms)
        self.__fdmt_plan = None

    @property
    def nuvrest(self):
        return self.fdmt_plan.nuvtotal // self.nuvwide

    @property
    def uv_shape(self):
        return (self.nuvrest, self.nt, self.ncin, self.nuvwide)


    @property
    def fdmt_plan(self):
        '''
        Lazy evanluate FDMT plan
        '''
        if self.__fdmt_plan is not None:
            return self.__fdmt_plan

        uvcells = self.uvcells
        nuvwide = self.nuvwide
        ncin = self.ncin

        self.__fdmt_plan = FdmtPlan(uvcells, self, prev_pipeline_plan=self.prev_plan)

        self.save_fdmt_plan_lut()
        
        if self.fdmt_plan.nuvtotal >= self.values.nuvmax:
            raise ValueError("Too many UVCELLS")

        self.upper_instructions = calc_grid_luts(self, True)
        self.lower_instructions = calc_grid_luts(self, False)
        self.save_grid_instructions(self.upper_instructions, 'upper')
        self.save_grid_instructions(self.lower_instructions, 'lower')
        self.upper_idxs, self.upper_shifts, self.lower_idxs, self.lower_shifts = calc_pad_lut(self, self.fft_ssr)
        self.save_pad_lut(self.upper_idxs, self.upper_shifts, 'upper')
        self.save_pad_lut(self.lower_idxs, self.lower_shifts, 'lower')

        # calculate DD grid id1 and offset values, there
        self.ddreader_config = make_ddreader_configs(self)

        # concatenate bytes of DMS and ddreader config
        self.ddreader_lut = np.frombuffer(self.dms.tobytes() + self.ddreader_config.tobytes(), dtype=np.uint32)
        self.save_lut(self.ddreader_lut, 'ddreader', 'value')


        return self.__fdmt_plan

        
    def save_lut(self, data, lutname, header, fmt='%d'):
        if self.values.uv is not None:
            filename = '{uvfile}.{lutname}.txt'.format(uvfile=self.values.uv, lutname=lutname)
            log.info('Saving {lutname} shape={d.shape} type={d.dtype} to {filename} header={header}'.format(lutname=lutname, d=data, filename=filename, header=header))
            np.savetxt(filename, data, fmt=fmt, header=header)

    def save_fdmt_plan_lut(self):
        fruns = self.fdmt_plan.runs
        d = []
        for irun, run in enumerate(fruns):
            if run is None:
                continue
            
            for icell, cell in enumerate(run.cells):
                if cell is None:
                    continue
                
                d.append([cell.a1,
                          cell.a2,
                          cell.uvpix[0],
                          cell.uvpix[1],
                          cell.chan_start,
                          cell.chan_end,
                          irun,
                          icell,
                          run.total_overlap,
                          run.max_idm,
                          run.max_offset,
                          run.offset_cff,
                          run.idm_cff,
                          run.fch1])

        d = np.array(d)

        header='ant1, ant2, u(pix), v(pix), chan1, chan2, irun, icell, total_overlap, max_idm, max_offset, offset_cff, idm_cff, fch1'
        fmt = '%d ' * 8 + ' %d '*3 + ' %f '*3
        self.save_lut(d, 'uvgrid.split', header, fmt=fmt)

        
    def save_grid_instructions(self, instructions, name):
        log.info('Got %d %s grid instructions', len(instructions), name)
        d = np.array([[i.target_slot, i.uvidx, i.shift_flag, i.uvpix[0], i.uvpix[1]] for i in instructions], dtype=np.int32)
        header ='target_slot, uvidx, shift_flag, upix, vpix'
        self.save_lut(d, 'gridlut.'+name, header)

    def save_pad_lut(self, idxs, shifts, name):
        d = np.array(idxs, dtype=np.int32)
        header ='upix, vpix, regidx'
        self.save_lut(d, 'padlut.'+name, header)

        d = np.array(shifts, dtype=np.int32)
        header = 'doshift'
        self.save_lut(d, 'doshift.'+name, header)

    def get_uv(self, blid):
        '''
        Returns the UV coordinates (in ns) of the given baseline ID
        
        :blid: Baseline ID
        :returns: (u, v) tuple in ns
        '''
        bldata = self.baselines(blid)
        return (bldata['UU'], bldata['VV'])

    @property
    def nf(self):
        '''Returns number of frequency channels'''
        return len(self.freqs)


    @property
    def fmin(self):
        '''
        Returns maximum frequency
        '''
        return self.freqs[0]

    @property
    def fmax(self):
        '''
        Returns minimum frequency
        '''
        return self.freqs[-1]

    @property
    def dmax(self):
        '''
        Returns maximum DM - placeholder for when we do DM gaps
        '''
        return max(self._dms)

    @property
    def nd(self):
        '''
        Returns numebr of dms
        '''
        return self._ndm

    @property
    def dms(self):
        return self._dms

    def __set_dms(self, in_dms):
        # TOOD: Check every DM is less than the maximum DM
        assert len(in_dms) <= self.values.max_ndm, f'Cant have more than {self.values.max_ndm} dms'
        self._ndm = len(in_dms)
        
        self._dms = np.zeros(self.values.max_ndm, dtype=np.uint32)
        self._dms[:len(in_dms)] = in_dms
        
        
    @property
    def tsamp_s(self):
        '''
        Returns sample time in seconds
        '''

        return self.__tsamp

    @property
    def nant(self):
        '''
        Returns number of antennas - basically from number of baselines and solving the quadratic equation
        '''
        na = int((1 + np.sqrt(1 + 8*self.nbl))/2)
        return na

    @property
    def maxant(self):
        '''
        Returns the largest antenna ID (1 based) in the baselines
        '''
        mx = max([max(craco.bl2ant(blid)) for blid in self.baseline_order])
        
        return mx


    def pointsource(self, coord, amp=1, noiseamp=0):
        '''
        Returns points source visibilities at given astropy skycoord
        :coord: Coordinate of point source
        :amp: amplitude
        :noiseamp: noise amplitude
        '''
        #psi = Angle('0.6d') # offset degrees - RA direction
        #theta = Angle('0.7d') # offset degrees - dec direction
        #expected_dec = plan.phase_center.dec + theta # not sure why I need a negative here
        # RA needs to be decremented by source cos dec
        #expected_ra = np.degrees(plan.phase_center.ra.rad - psi.rad/np.cos(expected_dec.rad))
        #expected_pos = SkyCoord(expected_ra, expected_dec, unit='deg')

        psi_rad= -coord.ra.rad*np.cos(coord.dec.rad) + plan.phase_center.rad
        
        lm = np.sin([psi.rad, theta.rad])

        ps = craco.pointsource(amp, lm, self.freqs, self.baseline_order, self.baselines, noiseamp)
        return ps
        
def add_arguments(parser):
    '''
    Add planning arguments
    '''
    parser.add_argument('--uv', help='Load antenna UVW coordinates from this UV file', default='uv_data')
    parser.add_argument('--pickle_fname', default='pipeline.pickle', help='File to dump and load pickle file')
    parser.add_argument('--npix', help='Number of pixels in image', type=int, default=256)
    parser.add_argument('--fov', help='Target field of view (angle with unit) - optionally comma separated for 2d. e.g. 0.9d or 1.1deg,0.9deg', default='1.1d')
    parser.add_argument('--os', help='Minimum number of pixels per beam. Optionaly comma separated for 2D e.g. 2.1 or 2.1,2.1', default='2.1')
    #parser.add_argument('--cell', help='Image cell size (arcsec). Overrides --os - possibly doesnt work')
    parser.add_argument('--nt', help='Number of times per block', type=int, default=256)
    parser.add_argument('--ndm', help='Number of DM trials', type=int, default=2)
    parser.add_argument('--max-ndm', help='Maximum number of DM trials. MUST AGREE WITH FIRMWARE - DO NOT CHANGE UNLESS YOU KNW WHAT YOUR DOING', type=int, default=1024)
    parser.add_argument('--max-nbl', help='Maximum number of baselines - cuts number of baselines to this value', type=int, default=36*35/2)
    parser.add_argument('--nbox', help='Number of boxcar trials', type=int, default=8)
    parser.add_argument('--boxcar-weight', help='Boxcar weighting type', choices=('sum','avg','sqrt'), default='sum')
    parser.add_argument('--nuvwide', help='Number of UV processed in parallel', type=int, default=8)
    parser.add_argument('--nuvmax', help='Maximum number of UV allowed.', type=int, default=8192-8+8) # For some reason NUREST is 1023 in craco_pybind11 but setting this to 8192 - 8 makes it hang - so put it back for now and be very very bloody careful.
    parser.add_argument('--ncin', help='Numer of channels for sub fdmt', type=int, default=32)
    parser.add_argument('--ndout', help='Number of DM for sub fdmt', type=int, default=186)
    parser.add_argument('--fdmt_scale', type=float, help='Scale FDMT output by this amount', default=1.0)
    parser.add_argument('--fft_scale', type=float, help='Scale FFT output by this amount. If both scales are 1, the output equals the value of frb_amp for crauvfrbsim.py', default=10.0)
    parser.add_argument('--show-image', action='store_true', help='Show image plots', default=False)
    parser.add_argument('--show-fdmt', action='store_true', help='Show FDMT plots', default=False)
    parser.add_argument('--save', action='store_true',  help='Save data as .npy for input, FDMT and image pipeline')
    parser.add_argument('--flag-ants', type=strrange, help='Ignore these 1-based antenna number')
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s', '--show', action='store_true', help='Show plots')
    parser.add_argument('--target-input-rms', type=float, default=512, help='Target input RMS')
    parser.add_argument('--calibration', bhelp='Calibration .bin file or root of Miriad files to apply calibration')


    # raelly only needed for pipeline
    parser.add_argument('--dflag-fradius', help='Dynamic flagging frequency radius. >0 to enable flagging', default=0, type=float)
    parser.add_argument('--dflag-tradius', help='Dynamic flagging time radius. >0 to enable flagging', default=0, type=float)
    parser.add_argument('--dflag-threshold', help='Dynamic flagging threshold. >0 to enable flagging', default=0, type=float)
    parser.add_argument('--dflag-tblk', help='Dynamif flagging block size. Must divide evenly into the block size (256 usually)', default=None, type=int)


def get_parser():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Plans a CRACO scan', formatter_class=ArgumentDefaultsHelpFormatter)
    add_arguments(parser)
    parser.set_defaults(verbose=False)
    return parser
               

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Plans a CRACO scan', formatter_class=ArgumentDefaultsHelpFormatter)
    add_arguments(parser)
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    log.info('Loading UV coordinates from file %s ', values.uv)
    f = uvfits.open(values.uv)
    plan = PipelinePlan(f, values)

    if values.show:
        f.plot_baselines()
        pylab.figure()
        nchans = [len(bc.freqs) for bc in plan.uvcells]
        pylab.hist(nchans, bins=50)
        pylab.xlabel('Number of channels per cell')
        pylab.ylabel('Number of cells')

        pylab.show()

    #dump_plan(plan, values.pickle_fname) # Large plans don't work - need to not serialise the bad bits
    #print((load_plan(values.pickle_fname)))
    
if __name__ == '__main__':
    _main()

    

    
