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
import craco
import uvfits
import warnings
import craco_kernels
from craco import triangular_index, make_upper
import fdmt

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


def get_uvcells(baselines, uvcell, freqs, Npix, plot=True):
    uvcells = []

    ucell, vcell = uvcell

    if plot:
        grid = np.zeros((Npix, Npix))

    # Updated for python 3
    #for blid, bldata in baselines.iteritems():
    for blid, bldata in baselines.items():
        #UU, VV WW are in seconds
        ulam = bldata['UU'] * freqs
        vlam = bldata['VV'] * freqs
        pix_offset = float(Npix)/2.0
        
        upix = np.round(ulam/ucell + pix_offset).astype(int)
        vpix = np.round(vlam/vcell + pix_offset).astype(int)
        if np.any((upix < 0) | (upix >= Npix) | (vpix < 0) | (vpix >= Npix)):
            warnings.warn('Pixel coordinates out of range')

        if plot:
            #pylab.plot(ulam/1e3, vlam/1e3)
            pylab.plot(ulam/ucell + pix_offset, vlam/vcell+pix_offset)

        uvpos = list(zip(upix, vpix))
        for istart, iend in craco.runidxs(uvpos):
            uvpix = uvpos[istart]
            assert uvpos[istart] == uvpos[iend]
            b = craco.BaselineCell(blid, uvpix, istart, iend, freqs[istart:iend+1], Npix)
            uvcells.append(b)
            grid[uvpix[1],uvpix[0]] += len(b.freqs)


    if plot:
        pylab.imshow(grid)
        pylab.show()

    uvcells = sorted(uvcells, key=lambda b:b.upper_idx)
    return uvcells

def calc_overlap_channels(chan_start, chan_end, minchan, ncin):
    '''
    Calculate overlap between 2 sets of channels:

    :chan_start: first channel of baseline. Must be >= minchan
    :chan_end: last_cahnenl of baseline (inclusive)
    :minchan: First channel in range
    :nchan: Number of chanenls in range

    >>> calc_overlap_channels(0,31,0,32)
    32

    >>> calc_overlap_channels(0,30,0,32)
    31

    >>> calc_overlap_channels(0,33,0,32)
    32

    >>> calc_overlap_channels(0+1,31+1,0,32)
    31

    >>> calc_overlap_channels(31,64,0,32)
    1

    >>> calc_overlap_channels(32,64,0,32)
    0

    >>> calc_overlap_channels(33,64,0,32)
    -1

    '''
    assert chan_end >= chan_start
    assert chan_start >= minchan, 'invalid minimum chan {} {}'.format(chan_start, minchan)
    maxchan = minchan + ncin # top end of channels - inclusive

    if chan_end < maxchan:
        overlap = chan_end - chan_start + 1
    else:
        overlap = maxchan - chan_start

    #print(chan_start, chan_end, minchan, maxchan, ncin, overlap)

    return overlap

def calc_overlap(blcell, minchan, ncin):
    return calc_overlap_channels(blcell.chan_start, blcell.chan_end, minchan, ncin)

def split_cells(cells, minchan, ncin):
    '''
    Given a list of UV cells and a channel range given by minchan and ncin, produces 2 lists
    The first list contains Cells that fint entirely within the chanel range. The second list
    contains cells outside the channel range

    :cells: List of UV cells
    :minchan: first channel in channel range
    :ncin: number of channels in range
    :returns: Tuple of lists
    '''
    full_cells = []
    leftover_cells = []
    endchan = minchan + ncin - 1 # final channel, inclusive

    for c in cells:
        assert c.chan_start >= minchan
        if c.chan_start >= minchan and c.chan_end <= endchan:
            full_cells.append(c)
        else:
            # split the cell in two parts
            npix = c.npix
            logging.debug('c %d-%d freq=%f-%f endchan=%d', c.chan_start, c.chan_end, c.freqs[0], c.freqs[-1], endchan)
            endidx = endchan - c.chan_start
            assert 0 <= endidx < ncin
            c1 = craco.BaselineCell(c.blid, c.uvpix, c.chan_start, endchan, c.freqs[:endidx+1], npix)
            c2 = craco.BaselineCell(c.blid, c.uvpix, endchan+1, c.chan_end, c.freqs[endidx+1:], npix)
            full_cells.append(c1)
            leftover_cells.append(c2)

    return full_cells, leftover_cells

class FdmtRun(object):
    def __init__(self, cells, plan):
        self.plan = plan
        mincell = min(cells, key=lambda cell:cell.chan_start)
        self.cells = cells
        self.chan_start = mincell.chan_start
        self.fch1 = mincell.fch1
        self.total_overlap = 0
        for uv in cells:
            overlap = calc_overlap(uv, self.chan_start, plan.pipeline_plan.ncin)
            logging.debug('Cell chan_start %s %s %s-%s overlap=%d', self.chan_start, uv, uv.chan_start, uv.chan_end, overlap)
            assert overlap > 0
            self.total_overlap += overlap

        assert self.max_idm <= plan.pipeline_plan.ndout, 'NDOUT is too small - needs to be at least %s' % self.max_idm

    @property
    def offset_cff(self):
        return craco_kernels.offset_cff(self.fch1, self.plan.pipeline_plan)

    @property
    def idm_cff(self):
        return craco_kernels.idm_cff(self.fch1, self.plan.pipeline_plan)

    @property
    def max_idm(self):
        dmax = self.plan.pipeline_plan.dmax
        return int(np.ceil(dmax*self.idm_cff))

    @property
    def max_offset(self):
        dmax = self.plan.pipeline_plan.dmax
        return int(np.ceil(dmax*self.offset_cff))

    def __str__(self):
        ncells = len(self.cells)
        return 'ncells={ncells} fch1={self.fch1} chan_start={self.chan_start} total_overlap={self.total_overlap}'.format(self=self, ncells=ncells)
        



class FdmtPlan(object):
    def __init__(self, uvcells, pipeline_plan):
        self.pipeline_plan = pipeline_plan
        nuvwide = self.pipeline_plan.nuvwide
        ncin = self.pipeline_plan.ncin
        uvcells_remaining = uvcells[:] # copy array
        fdmt_runs = []
        run_chan_starts = []
        run_fch1 = []
        runs = []
        while len(uvcells_remaining) > 0:
            logging.debug('Got %d/%d uvcells remaining', len(uvcells_remaining), len(uvcells))
            minchan = min(uvcells_remaining, key=lambda uv:(uv.chan_start, uv.blid)).chan_start
            possible_cells = filter(lambda uv:calc_overlap(uv, minchan, ncin) > 0, uvcells_remaining)
            #possible_cells_copy = []

            #possible_cells_copy[:] = possible_cells
            
            #logging.debug('Got %d possible cells', len(possible_cells))

            # sort as best we can so that it's stable - I.e. we get hte same answer every time
            best_cells = sorted(possible_cells, key=lambda uv:(calc_overlap(uv, minchan, ncin), uv.blid, uv.upper_idx), reverse=True)
            logging.debug('Got %d best cells. Best=%s overlap=%s', len(best_cells), best_cells[0], calc_overlap(best_cells[0], minchan, ncin))
            used_cells = best_cells[0:min(nuvwide, len(best_cells))]
            full_cells, leftover_cells = split_cells(used_cells, minchan, ncin)
            run = FdmtRun(full_cells, self)
            run_chan_starts.append(run.chan_start)
            run_fch1.append(run.fch1)
            fdmt_runs.append(full_cells)
            runs.append(run)
            total_overlap = run.total_overlap
            #logging.debug('minchan=%d npossible=%d used=%d full=%d leftover=%d total_overlap=%d', minchan, len(possible_cells), len(used_cells), len(full_cells), len(leftover_cells), total_overlap)
            
            # Remove used cells
            uvcells_remaining = [cell for cell in uvcells_remaining if cell not in used_cells]

            # Add split cells
            uvcells_remaining.extend(leftover_cells)
            
        nruns = len(fdmt_runs)
        nuvtotal = nruns*nuvwide

        ndout = self.pipeline_plan.ndout
        nd = self.pipeline_plan.nd
        nt = self.pipeline_plan.nt
        #square_history_size = ndout*nuvtotal*(nt + nd)
        square_history_size = sum(nuvwide*(nd + nt)*ndout for run in runs)
        minimal_history_size = sum(nuvwide*(run.max_offset+ nt)*run.max_idm for run in runs)
        efficiency = float(len(uvcells))/float(nuvtotal)
        required_efficiency = float(nuvtotal)/8192.0
        
        logging.info('FDMT plan has ntotal=%d of %d runs with packing efficiency %f. Grid read requires efficiency of > %f of NUV=8192. History size square=%d minimal=%d =%d 256MB HBM banks', nuvtotal, nruns, efficiency, required_efficiency, square_history_size, minimal_history_size, minimal_history_size*4/256/1024/1024)
        self.fdmt_runs = fdmt_runs
        self.run_chan_starts = run_chan_starts
        self.run_fch1 = run_fch1
        self.runs = runs

        # create an FDMT object for each run so  we can use it to calculate the lookup tbales
        #     def __init__(self, f_min, f_off, n_f, max_dt, n_t, history_dtype=None):

        fdmts = [fdmt.Fdmt(fch1, self.pipeline_plan.foff, ncin, self.pipeline_plan.ndout, 1) for fch1 in self.run_fch1]
        fdmt_luts = np.array([thefdmt.calc_lookup_table() for thefdmt in fdmts])
        niter = int(np.log2(ncin))
        assert 1<<niter == ncin
        # final LUTs we need to copy teh same LUT for every NUVWIDE
        assert fdmt_luts.shape == (nruns, niter, ncin-1, 2)
        self.fdmt_lut = np.repeat(fdmt_luts[:,:,:,np.newaxis, :], nuvwide, axis=3)
        expected_lut_shape = (nruns, niter, ncin-1, nuvwide, 2)
        assert self.fdmt_lut.shape == expected_lut_shape, 'Incorrect shape for LUT=%s expected %s' % (self.fdmt_lut.shape, expected_lut_shape)
        

        self.nruns = nruns
        self.nuvtotal = nuvtotal
        self.total_nuvcells = sum([len(p) for p in fdmt_runs])

        # find a cell with zero in it
        self.zero_cell = None
        for irun, run in enumerate(self.runs):
            if len(run.cells) < nuvwide:
                self.zero_cell = (irun, len(run.cells))

        if self.zero_cell is None:
            self.zero_cell = (len(self.runs), 0)

        assert self.zero_cell != None
        assert self.zero_cell[0] < self.pipeline_plan.nuvrest, 'Not enough room for FDMT zero cell'
        assert self.zero_cell[1] < nuvwide

        assert self.zero_cell[0] < self.nruns
        assert self.zero_cell[1] < ncin
        #assert self.get_cell(self.zero_cell) != None

        logging.info("FDMT zero cell is %s=%s", self.zero_cell, self.zero_cell[0]*nuvwide+self.zero_cell[1])

        uvmap = {}
        for irun, run in enumerate(self.runs):
            for icell, cell in enumerate(run.cells):
                irc = (irun, icell)
                uvmap[cell.uvpix_upper] = uvmap.get(cell.uvpix_upper, [])
                uvmap[cell.uvpix_upper].append(irc)

        self.__uvmap = uvmap

    def cell_iter(self):
        '''
        Iteration over all the cells
        '''
        for run in self.runs:
            for cell in run.cells:
                yield cell

    def find_uv(self, uvpix):
        '''
        Returns the run and cell index of FDMT Cells that have the given UV pixel

        uvpix must be upper hermetian
        
        Returns a list of tuples
        irun = run index
        icell = cellindex inside the run

        You can find the cell with 
        self.get_cell((irun, icell))
        
        :uvpix: 2-tuple of (u, v)
        :returns: 2-typle (irun, icell)
        '''

        assert uvpix[0] >= uvpix[1], 'Uvpix must be upper hermetian'
        # speed optimisation
        cell_coords2 = self.__uvmap.get(uvpix,[])
        
        #cell_coords = []
        #for irun, run in enumerate(self.runs):
        #    for icell, cell in enumerate(run.cells):
        #        if cell.uvpix_upper == uvpix:
        #            cell_coords.append((irun, icell))


        #print(uvpix, 'Version1', cell_coords, 'Verion2', cell_coords2)
        #assert cell_coords == cell_coords2

        return cell_coords2

    def get_cell(self, cell_coord):
        irun, icell = cell_coord
        if cell_coord == self.zero_cell:
            return None
        else:
            return self.runs[irun].cells[icell]

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
        sorter = lambda c: triangular_index(c[0],c[1], plan.npix)
    else:
        uvpix_list = [cell.uvpix_lower for cell in uvcells if cell.uvpix[0] != cell.uvpix[1]] # don't include diagonal
        sorter = lambda c: triangular_index(c[1],c[0], plan.npix, raster='yx') # triangular_index requires upper_triangular coordinates
        
    unique_uvcoords = set(uvpix_list)
    # sorts in raster order

    unique_uvcoords = sorted(unique_uvcoords, key=sorter)
    ncoord = len(unique_uvcoords)
    logging.info('Got %d unique UV coords. Upper=%s', ncoord, upper)

    # check
    check = False
    if check:
        grid = np.zeros((plan.npix, plan.npix))
        for i, (u, v) in enumerate(unique_uvcoords):
            grid[v,u] = i+1

        pylab.imshow(np.ma.masked_equal(grid, 0), interpolation='none')
        pylab.show()
    
    fplan = plan.fdmt_plan
    fruns = fplan.runs
    remaining_fdmt_cells = [] # this is ana rray we keep for self-testing to make sure we used everything
    for run in fruns:
        remaining_fdmt_cells.extend(run.cells)

    ngridreg = plan.ngridreg
    all_instructions = []

    nwrites = (ncoord + ngridreg - 1) // ngridreg # number of writes
    assert nwrites <= 4096*2, 'Not enough clocks to write data in! nwrites={}'.format(nwrites)
    logging.info('Need to write %d groups of %d register to pad function', nwrites, ngridreg)

    #for iwrite in xrange(nwrites):
    for iwrite in range(nwrites):
        # Add available cells
        n = min(ngridreg, ncoord - iwrite*ngridreg)

        slot_uv_coords = unique_uvcoords[iwrite*ngridreg:iwrite*ngridreg + n]
        instructions = []

        for islot, slot_uv_coord in enumerate(slot_uv_coords):
            logging.debug('Considering islot=%d slot=%s', islot, slot_uv_coord)
            # The FDMT has all its data in the upper hermetian - so we always use upper hermetian coordinates
            upper_slot_uv_coord = make_upper(slot_uv_coord)
             # Find irun and icell 
            for cell_coord in fplan.find_uv(upper_slot_uv_coord):
                # cell_coord is the location in the FDMT output buffer
                cell = fplan.get_cell(cell_coord) # this is the actual baseline cell. Don't use it other than to remove it form our check buffer
                logging.debug('removing cell %s %s', cell, cell_coord)
                remaining_fdmt_cells.remove(cell)

                inst = AddInstruction(plan, islot, cell_coord, slot_uv_coord)
                instructions.append(inst)

        # Grid reading reads 2 UV/clk, so we need to add a dummy add instruction
        if len(instructions) % 2 == 1: # if odd
            instructions.append(AddInstruction(plan, ngridreg-1, fplan.zero_cell, (-1,-1)))

        # See shift marker for last 2 instruction, as we read 2UV/clk and I don't know which mattes
        #instructions[-2].shift = True
        instructions[-1].shift = True

        all_instructions.extend(instructions)
        logging.debug('Added %d instructions with %d available', len(instructions), len(remaining_fdmt_cells))

    # Check instructions

    assert all_instructions[-1].shift == True
    #assert all_instructions[-2].shift == True
    if upper:
        assert len(remaining_fdmt_cells) == 0

    num_shifts =  sum(map(lambda inst:inst.shift == True, all_instructions))

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
    upper_inst = plan.upper_instructions
    lower_inst  = plan.lower_instructions
    upper_inputs = get_pad_input_registers(upper_inst)
    lower_inputs = get_pad_input_registers(lower_inst)
    uvpix_set = set([cell.uvpix_upper for cell in plan.uvcells])

    lut = []
    upper_registers = []
    lower_registers = []

    # Update for python 3
    #upper_registers.extend(upper_inputs.next())
    #upper_registers.extend(upper_inputs.next())
    #
    #lower_registers.extend(lower_inputs.next())
    #lower_registers.extend(lower_inputs.next())

    upper_registers.extend(next(upper_inputs))
    upper_registers.extend(next(upper_inputs))
                                            
    lower_registers.extend(next(lower_inputs))
    lower_registers.extend(next(lower_inputs))
    
    upper_shifts = []
    lower_shifts = []
    upper_idxs = []
    lower_idxs = []

    npix = plan.npix
    
    #for v in xrange(npix):
    for v in range(npix):
        #for ublk in xrange(npix/ssr):
        for ublk in range(npix//ssr):
            lower_shift = False
            upper_shift = False
            #for iu in xrange(ssr):
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
                        i = -1 # which means 0

                    upper_idxs.append((u,v,i))
                else: # lower hermeitan
                    if (v,u) in uvpix_set:
                        i = lower_registers.index(uv)
                        lower_shift |= i >= ssr - 1
                    else:
                        i = -1
                        
                    lower_idxs.append((u,v,i))


            if upper_shift:
                upper_registers[:ssr] = upper_registers[ssr:]
                try:
                    # update for python 3
                    #upper_registers[ssr:] = upper_inputs.next()
                    upper_registers[ssr:] = next(upper_inputs)
                except StopIteration:
                    upper_finished = True


            if lower_shift:
                lower_registers[:ssr] = lower_registers[ssr:]
                try:
                    # Update for python 3
                    #lower_registers[ssr:] = lower_inputs.next()
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
    def __init__(self, f, values):
        self.values = values
        logging.info('making Plan values=%s', values)

        umax, vmax = f.get_max_uv()
        lres, mres = 1./umax, 1./vmax
        baselines = f.baselines
        nbl = len(baselines)
        freqs = f.channel_frequencies

        # Cant handle inverted bands - this is assumed all over the code. It's boring
        assert freqs.min() == freqs[0]
        assert freqs.max() == freqs[-1]
        Npix = values.npix

        if values.cell is not None:
            lcell, mcell = map(craco.arcsec2rad, values.cell.split(','))
            los, mos = lres/lcell, mres/mcell
        else:
            los, mos = map(float, values.os.split(','))
            lcell = lres/los
            mcell = mres/mos
            
        lfov = lcell*Npix
        mfov = mcell*Npix
        ucell, vcell = 1./lfov, 1./mfov
        fmax = freqs.max()
        foff = freqs[1] - freqs[0]
        lambdamin = 3e8/fmax
        umax_km = umax*lambdamin/1e3
        vmax_km = vmax*lambdamin/1e3
        
        logging.info('Nbl=%d Fch1=%f foff=%f nchan=%d lambdamin=%f uvmax=%s max baseline=%s resolution=%sarcsec uvcell=%s arcsec uvcell= %s lambda FoV=%s deg oversampled=%s',
                 nbl, freqs[0], foff, len(freqs), lambdamin, (umax, vmax), (umax_km, vmax_km), np.degrees([lres, mres])*3600, np.degrees([lcell, mcell])*3600., (ucell, vcell), np.degrees([lfov, mfov]), (los, mos))

        #print(baselines)
        
        uvcells = get_uvcells(baselines, (ucell, vcell), freqs, Npix)
        logging.info('Got Ncells=%d uvcells', len(uvcells))
        d = np.array([(v.a1, v.a2, v.uvpix[0], v.uvpix[1], v.chan_start, v.chan_end) for v in uvcells], dtype=np.int32)
        np.savetxt(values.uv+'.uvgrid.txt', d, fmt='%d',  header='ant1, ant2, u(pix), v(pix), chan1, chan2')

        self.uvcells = uvcells
        self.nd = values.ndm
        self.nt = values.nt
        self.freqs = freqs
        self.npix = Npix
        self.nbox = values.nbox
        self.boxcar_weight = values.boxcar_weight
        self.nuvwide = values.nuvwide
        self.nuvmax = values.nuvmax
        assert self.nuvmax % self.nuvwide == 0
        self.nuvrest = self.nuvmax // self.nuvwide
        self.ncin = values.ncin
        self.ndout = values.ndout
        self.foff = foff
        self.dtype = np.complex64 # working data type
        self.threshold = values.threshold
        self.nbl = nbl
        self.fdmt_scale = self.values.fdmt_scale
        self.fft_scale = self.values.fft_scale
        self.fft_ssr = 16 # number of FFT pixels per clock - "super sample rate"
        self.ngridreg = 16 # number of grid registers to do
        assert self.threshold >= 0, 'Invalid threshold'
        self.fdmt_plan = FdmtPlan(uvcells, self)
        self.save_fdmt_plan_lut()

        
        if self.fdmt_plan.nuvtotal >= values.nuvmax:
            raise ValueError("Too many UVCELLS")

        self.upper_instructions = calc_grid_luts(self, True)
        self.lower_instructions = calc_grid_luts(self, False)
        self.save_grid_instructions(self.upper_instructions, 'upper')
        self.save_grid_instructions(self.lower_instructions, 'lower')
        self.upper_idxs, self.upper_shifts, self.lower_idxs, self.lower_shifts = calc_pad_lut(self, self.fft_ssr)
        self.save_pad_lut(self.upper_idxs, self.upper_shifts, 'upper')
        self.save_pad_lut(self.lower_idxs, self.lower_shifts, 'lower')

        filehandler = open("pipeline.obj", 'wb') 
        pickle.dump(self, filehandler)
        filehandler.close()
        
    def save_lut(self, data, lutname, header, fmt='%d'):
        filename = '{uvfile}.{lutname}.txt'.format(uvfile=self.values.uv, lutname=lutname)
        logging.info('Saving {lutname} shape={d.shape} type={d.dtype} to {filename} header={header}'.format(lutname=lutname, d=data, filename=filename, header=header))
        np.savetxt(filename, data, fmt=fmt, header=header)

    def save_fdmt_plan_lut(self):
        fruns = self.fdmt_plan.runs
        d = []
        for irun, run in enumerate(fruns):
            for icell, cell in enumerate(run.cells):
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
        logging.info('Got %d %s grid instructions', len(instructions), name)
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
        return self.nd
            
            
        

def add_arguments(parser):
    '''
    Add planning arguments
    '''
    parser.add_argument('--uv', help='Load antenna UVW coordinates from this UV file')
    parser.add_argument('--npix', help='Number of pixels in image', type=int, default=256)
    parser.add_argument('--os', help='Number of pixels per beam', default='2.1,2.1')
    parser.add_argument('--cell', help='Image cell size (arcsec). Overrides --os')
    parser.add_argument('--nt', help='Number of times per block', type=int, default=16)
    parser.add_argument('--ndm', help='Number of DM trials', type=int, default=16)
    parser.add_argument('--nbox', help='Number of boxcar trials', type=int, default=4)
    parser.add_argument('--boxcar-weight', help='Boxcar weighting type', choices=('sum','avg','sqrt'), default='sum')
    parser.add_argument('--nuvwide', help='Number of UV processed in parallel', type=int, default=8)
    parser.add_argument('--nuvmax', help='Maximum number of UV allowed.', type=int, default=8192)
    parser.add_argument('--ncin', help='Numer of channels for sub fdmt', type=int, default=32)
    parser.add_argument('--ndout', help='Number of DM for sub fdmt', type=int, default=32)
    parser.add_argument('--threshold', type=float, help='Threshold for candidate grouper', default=10)
    parser.add_argument('--fdmt_scale', type=float, help='Scale FDMT output by this amount', default=1.0)
    parser.add_argument('--fft_scale', type=float, help='Scale FFT output by this amount. If both scales are 1, the output equals the value of frb_amp for crauvfrbsim.py', default=10.0)
    parser.add_argument('--show-image', action='store_true', help='Show image plots')
    parser.add_argument('--show-fdmt', action='store_true', help='Show FDMT plots')
    parser.add_argument('--save', action='store_true',  help='Save data as .npy for input, FDMT and image pipeline')


def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Plans a CRACO scan', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('-s', '--show', action='store_true', help='Show plots')
    add_arguments(parser)
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    logging.info('Loading UV coordinates from file %s ', values.uv)
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

    filename = "pipeline.obj"
    filehandler = open(filename, 'rb') 
    object = pickle.load(filehandler)
    print(object.values)

if __name__ == '__main__':
    _main()
