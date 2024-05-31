#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2020
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
from . import craco
from . import fdmt
from . import craco_kernels
from .craco import BaselineCell


log = logging.getLogger(__name__)

__author__ = "Keith Bannister <keith.bannister@csiro.au>"


def find_cells(uvcells, blid):
    '''
    Returns a list of cells with th given baseline ID
    Sorted by channel_start
    Returns empty list if blid not found which can happen for U=0 V=0 bug
    '''
    cells =[cell for cell in uvcells if cell is not None and cell.blid == blid] 
    cells.sort(key=lambda c:c.chan_start)
    
    return cells

def flat_cell_iter(runs):
    '''
    Flattens all the runs and returns a 1D generator of all non-None cells
    '''
    for run in runs:
        if run is not None:
            for cell in run.cells:
                if cell is not None:
                    yield cell
                        

def baselinestr(c):
    '''
    Print summary of a baselineCell
    '''
    if c is None:
        s = f'------ (--,--)'
    else:
        s = f'{c.chan_start}-{c.chan_end}={c.nchan} {c.uvpix_upper}'
    return s
                        
def rangestr(i, cells):
    
    if i >= len(cells):
        s = 'XXX-XXX (XXX,XXX)'
    else:
        c = cells[i]
        s = baselinestr(c) 
    return s

def channel_overlap(c1_start, c1_end, c2_start, c2_end):
    '''
    Returns the overlap of 2 intervals: c1_start:c1_end (inclusive) and c2_start:c2_end (inclusive)
    Returns 0 if disjoint
    '''
    cstart = max(c1_start, c2_start)
    cend = min(c1_end, c2_end)    
    overlap = max(cend - cstart + 1, 0)
    
    #print('overlap', c1, c2, cstart, cend, overlap)
    return overlap
    

def cell_channel_overlap(c1, c2):
    '''
    Returns the number of overlapping channels between cell 1 and cell 2
    '''
    return channel_overlap(c1.chan_start, c1.chan_end, c2.chan_start, c2.chan_end)


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
    maxchan = minchan + ncin # top end of channels - inclusivetoi

    if chan_end < maxchan:
        overlap = chan_end - chan_start + 1
    else:
        overlap = maxchan - chan_start

    #print(chan_start, chan_end, minchan, maxchan, ncin, overlap)

    return overlap

def calc_overlap(blcell, minchan, ncin):
    return calc_overlap_channels(blcell.chan_start, blcell.chan_end, minchan, ncin)
    
def get_cell_with_highest_overlap(cell, last_cells):
    return max(last_cells, key=lambda c: cell_channel_overlap(cell, c))

def chan_with_best_overlap(mincell, chan_positions, chan_width):
       # find channel range that has the largest overlap with this cell
    overlaps = np.array([calc_overlap_channels(mincell.chan_start, 
                                               mincell.chan_end, 
                                               c, 
                                               chan_width) if c <= mincell.chan_start <= c+chan_width-1 else 0 
                         for c in chan_positions])
    assert len(overlaps) == len(chan_positions)
    best_overlap_idx = np.argmax(overlaps)  
    run_start_chan = chan_positions[best_overlap_idx]
    return run_start_chan

def fixed_algorithm(cells, chan_positions, chan_width=32, nuvwide=8):
    '''
    Start with a fixed set of given by the chan_positions array.
    When starting a run, find an appropriate channel from the chan_positions array 
    and use that
    '''
    runs = []
    run_start_chans = []
        
    cells_left = cells[:] # copy list
    minimum_nruns = (len(cells) + 8 - 1) // 8 # // in python means integer division
    log.info(f'Started with {len(cells)} cells which would need at least {minimum_nruns} runs')      
    while len(cells_left) > 0:
        # Find start channel from list of chan_positions that has. We have to start somewhere.
        mincell = min(cells_left, key=lambda cell:cell.chan_start)        
        
        # find that start channel
        lowest_chan = mincell.chan_start        
        
        method = 2

        if method == 1:
            # make list of channels from chan_positions that brackets lowest_chan

            possible_chans = list([c for c in chan_positions if c <= lowest_chan <= c+chan_width-1])
            assert len(possible_chans) > 0, 'Channel positions cant handle channel {}'.format(lowest_chan)
        
            # channel to use is the largest of these - which gives the best chance of including lots of other channels
        
            run_start_chan = max(possible_chans)
        else:
            run_start_chan = chan_with_best_overlap(mincell, chan_positions, chan_width)
 
        
        assert run_start_chan in chan_positions            
        
        # end channel of this run is
        run_end_chan = run_start_chan + chan_width-1 # inclusive
        
        # make a list of all UVs with a positive overlap - these UVs could be in this run.
        possible_uvs = [cell for cell in cells_left 
                        if calc_overlap_channels(cell.chan_start, 
                                                 cell.chan_end, 
                                                 run_start_chan, 
                                                 chan_width) > 0]
        
        # Sort UVs by how much overlap they have with this channel range
        # UVs with the largest overlap will be first in this list
        uv_sorted = sorted(possible_uvs, key=lambda cell:calc_overlap_channels(cell.chan_start, cell.chan_end, run_start_chan, chan_width), reverse=True)               
        
        # Extract no more than 8 Uvs from the list of UVs with the most overlap
        this_run_uvs = uv_sorted[0:min(nuvwide, len(uv_sorted))]
        
        uvs_for_this_run = []
        
        
        # for each UV in this run
        for cell in this_run_uvs:
            
            # remove it from the list of things to do
            cells_left.remove(cell)            
                      
            # If this uv final channel is beynd the channel range in this run we have to split it in two pieces
            uv_start_chan = cell.chan_start
            uv_end_chan = cell.chan_end
            assert uv_end_chan >= uv_start_chan
            if uv_end_chan > run_end_chan:
                # Add the subset of channels we *can* use to this run
                cell1, cell2 = cell.split_at(run_end_chan)
                uvs_for_this_run.append(cell1)
                
                # Add the leftover channels to the list of things to do
                cells_left.append(cell2)                
                
            else: # The whole channel range does fit in the run
                # Add the whole UV to this run
                uvs_for_this_run.append(cell)
                
        runs.append((run_start_chan, uvs_for_this_run))
                

    increase = len(runs)/minimum_nruns
    log.info(f'With {len(chan_positions)} channels at {chan_positions} we needed {len(runs)} runs. An increase of {increase:0.3f} ')
            
    return runs

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
            log.debug('c %d-%d freq=%f-%f endchan=%d', c.chan_start, c.chan_end, c.freqs[0], c.freqs[-1], endchan)
            endidx = endchan - c.chan_start
            assert 0 <= endidx < ncin
            c1 = craco.BaselineCell(c.blid, c.uvpix, c.chan_start, endchan, c.freqs[:endidx+1], npix)
            c2 = craco.BaselineCell(c.blid, c.uvpix, endchan+1, c.chan_end, c.freqs[endidx+1:], npix)
            full_cells.append(c1)
            leftover_cells.append(c2)

    return full_cells, leftover_cells

def fixed_create_fdmt_runs(pipeline_plan, fdmt_plan, uvcells, nuvwide, ncin):
        ncin = pipeline_plan.ncin
        chan_positions = np.arange(0,pipeline_plan.nf+ncin,ncin) # just need to have more than enough channels
        run_cells = fixed_algorithm(pipeline_plan.uvcells, chan_positions, chan_width=32)
        runs = [FdmtRun.from_cells(cells, fdmt_plan, chan_start) for (chan_start, cells) in run_cells]

        return runs
        
def greedy_create_fdmt_runs(theplan, uvcells, nuvwide, ncin):
    uvcells_remaining = uvcells[:] # copy array
    runs = []
    while len(uvcells_remaining) > 0:
        log.debug('Got %d/%d uvcells remaining', len(uvcells_remaining), len(uvcells))
        minchan = min(uvcells_remaining, key=lambda uv:(uv.chan_start, uv.blid)).chan_start
        possible_cells = [uv for uv in uvcells_remaining if calc_overlap(uv, minchan, ncin) > 0]

        # Do not know how to get a length of iterator in python3, comment it out here
        #log.debug('Got %d possible cells', len(possible_cells))

        # sort as best we can so that it's stable - I.e. we get hte same answer every time
        best_cells = sorted(possible_cells, key=lambda uv:(calc_overlap(uv, minchan, ncin), uv.blid, uv.upper_idx), reverse=True)
        log.debug('Got %d best cells. Best=%s overlap=%s', len(best_cells), best_cells[0], calc_overlap(best_cells[0], minchan, ncin))
        used_cells = best_cells[0:min(nuvwide, len(best_cells))]
        full_cells, leftover_cells = split_cells(used_cells, minchan, ncin)
        run = FdmtRun.from_cells(full_cells, theplan)
        runs.append(run)
        # create lookup table for each run
            
        total_overlap = run.total_overlap
        # Do not know how to get a length of iterator in python3, comment it out here
        #log.debug('minchan=%d npossible=%d used=%d full=%d leftover=%d total_overlap=%d', minchan, len(possible_cells), len(used_cells), len(full_cells), len(leftover_cells), total_overlap)
            
        # Remove used cells
        uvcells_remaining = [cell for cell in uvcells_remaining if cell not in used_cells]

        # Add split cells
        uvcells_remaining.extend(leftover_cells)

    return runs

def make_fdmt_lut(freqs, foff, ncin, ndout, nuvwide, do_correction=True):
    '''
    Make FDMT lookup tables suitable for hardware for given list of frequencies for each run
    :freqs: center chanenlf frequencies
    :foff:channel interval
    :ncin: number of channels per FDMT
    :ndout: output number of dm trials
    :nuvwide: number of FDMTs per run
    :do_correction: true if you want to the the clever FDMT correaction (you want True)
    '''
    
    # fill in a list of channels from the runs
    # if the run is none just use the fch1 from the pipeline
    # The FDMT will execute but zeros will be in there anyway
    # Those FDMTs will be executed by they'll have zero
    # do_correction = True makes the DM track pretty thign
    # do_correction = False basically doubles the width which makes the whole thing very wide and noisey
    # SEe "Testing image pipeline with impulses.ipynb.

    #freqs = [self.pipeline_plan.fmin if run is None else run.fch1 for run in runs]
    fdmts = [fdmt.Fdmt(freq-foff/2.0, foff, ncin, ndout, 1, do_correction=do_correction) for freq in freqs]
    fdmt_luts = np.array([thefdmt.calc_lookup_table() for thefdmt in fdmts])
    niter = int(np.log2(ncin))
    # final LUTs we need to copy teh same LUT for every NUVWIDE
    nruns = len(freqs)
    assert fdmt_luts.shape == (nruns, ncin-1, 2)
    fdmt_lut = np.repeat(fdmt_luts[:,:,np.newaxis, :], nuvwide, axis=2)
    expected_lut_shape = (nruns, ncin-1, nuvwide, 2)
    assert fdmt_lut.shape == expected_lut_shape, f'Incorrect shape for LUT={self.fdmt_lut.shape}  expected {expected_lut_shape}'

    return fdmt_lut

class FdmtPlanContainer:
    '''
    This is a utility class to help organise my brain. it matains a list of plan.NUVMAX FdmtRuns.
    If the run is active the list contains the run, otherwise it contains None
    If the container is initialised from a previous container, it sets up empty runs with the same chan_start
    at all the previous locations and expects to be poked and prodded until it's full
    I'm not entirely happy with the design of this (In particular, i painted myself into a corner to working out 
    how to create a new Run - it needs extra parameters) - so it isn't generic. It relised on the fixed algorithm
    way of allocating runs. I'm sure it can be improved.
    
    '''
    def __init__(self, plan:'FdmtPlan', runs:list, prev_plan:'FdmtPlanContainer'):
        '''
        Creates the container. If prev_plan is None, it uses the list of FdmtRuns to intiialise the internal list
        If prev_plan is not None it uses the runs from the previous plan but it creates new ones with the same bottom channel
        but no baseline cells, so they can be added later.
        '''
        self.nuvmax = plan.pipeline_plan.nuvrest_max
        self.nuvwide = plan.nuvwide
        self.__runs = [None for r in range(self.nuvmax)]
        self.plan = plan
        self.uvcells = list(flat_cell_iter(runs))
        self.prev_plan = prev_plan
        self.total_overlap = 0
        self.chan_positions = np.arange(0,32*12,32)
        self.chan_width = 32
        
        if prev_plan is None:
            self.set_runs(runs) # copy runs and cells configured
        else:
            for irun, run in enumerate(prev_plan.__runs):
                if run is not None:
                    self.set_run(irun, run.chan_start) # make empty run for every cell where appropriate          
                    
    def set_runs(self, runs):
        '''
        Initialise with the given runs
        '''
        for irun, run in enumerate(runs):
            self.__runs[irun] = run
            
        assert len(self.__runs) == self.nuvmax

    def cell_array_of(self, func):
        a = np.zeros((self.nuvmax, self.nuvwide))
        for u in range(self.nuvmax):
            for w in range(self.nuvwide):
                c = self.cell_at((u,w))
                n = np.nan if c is None else func(c)
                a[u,w] = n

        return a

    @property
    def nchan_array(self):
        return self.cell_array_of(lambda c:c.nchan)

    @property
    def start_chan_array(self):
        return np.array([np.nan if r is None else r.chan_start for r in self.__runs])

    @property
    def ncell_array(self):
        return np.array([np.nan if r is None else r.ncell for r in self.__runs])
                    
    @property
    def runs(self):
        '''
        Returns a copy of the list so people dont meddle with it
        '''
        return list(self.__runs)

    @property
    def valid_runs(self):
        '''
        Returns a list of runs but only up the last valid run index
        '''
        return self.__runs[:self.last_valid_run_index+1]
    
    @property
    def nruns(self):
        '''
        Returns the number of valid runs - will be <= nuvmax
        '''
        return sum(r is not None for r in self.__runs)
    
    def cell_iter(self):
        return flat_cell_iter(self.__runs)        
        
    @property
    def blids(self):
        return set([uv.blid for uv in self.uvcells])

    @property
    def nchan(self):
        '''
        Returns the total number of channels allocated in allcells
        '''
        return sum(c.nchan for c in self.cell_iter())
        
    def uvcells_of(self, blid):
        '''
        Return all uvcells in the runs with the given baseline id
        '''
        return find_cells(self.cell_iter(), blid)
        
    def run_at(self, irun:int) -> 'FdmtRun':
        '''
        Returns the run at the given index
        '''
        
        if irun >= len(self.__runs):
            raise IndexError(f'Invalid index {irun} of {len(self.__runs)}')
    
        return self.__runs[irun]
    
    def cell_at(self, coord) -> BaselineCell:
        '''
        Returns the BaselineCell at the given coordinate. 
        Coordinate is the (irun, icell) as a tuple
        If the run or cell are None it returns None
        '''
        irun, icell = coord
        
        run = self.run_at(irun)
        cell = None
        if run is not None:
            if icell >= len(run.cells):
                raise IndexError(f'Invalid cell index {icell} of {len(run.cells)} {run.cells}')
                
            cell = run.cells[icell]
                                                         
        return cell
    
    
    def set_run(self, irun, chan_start):
        '''
        Create a new run at the given index with the given chan_start
        If there's already a run at that index with the given chan_start it does nothing
        If the run has a different chan_start it throws an exception
        '''
        r = self.run_at(irun)
        if r is None:
            r = FdmtRun(self.plan, chan_start)
            self.__runs[irun] = r
            log.debug(f'Added run at {irun} with chan_start={chan_start} prev={None if self.prev_plan is None else self.prev_plan.run_at(irun)}')
        else:
            assert r.chan_start == chan_start, f'Inavlid chan start {r.chan_start} != {chan_start} {r}'
            
        return r
    
    @property
    def last_valid_run_index(self):
        '''
        Returns the index of the last valid run
        Equivalent to the number of UVWIDE we need to compute - 1
        '''
        idx = None
        for irun, run in enumerate(self.__runs):
            if run is not None:
                idx = irun
        
        assert idx is None or self.__runs[idx] is not None
        
        return idx
    
    @property
    def next_empty_run_index(self):
        '''
        Returns the index where a new run can be allocated.
        If there's a previous plan, then a run can be allocated where there is no current run and no previously planned run
        '''
        for i in range(len(self.__runs)):
            if self.run_at(i) is None and self.prev_plan is not None and self.prev_plan.run_at(i) is None:
                break
        return i
    
    def add_run(self, chan_start):
        '''
        Create new run at the first available slot
        '''

        idx = self.next_empty_run_index
        r = self.set_run(idx, chan_start)
        return r, idx
    
    def set_cell_coord(self, cell, coord):
        '''
        Insert a baselineCell at the given coordinate. Coordinate is (irun, icell)
        If there's  aprevious plan it will check if the previous cell had the same baseline ID
        and the new cell has an overlap of >0
        '''
        prev_cell = None
        overlap = 0
        if self.prev_plan is not None: # that it's sensible to assign this cell to this fdmt run/index
            prev_cell = self.prev_plan.cell_at(coord)
            if prev_cell is not None: # It's OK to run over an empty cell
                assert cell.blid == prev_cell.blid, 'Trying to set run for different BLID'
                overlap = cell_channel_overlap(cell, prev_cell)
                
                # There can sometimes be zero overlap, I guess we can't do much about that
                assert overlap > 0, f'Zero overlap for {coord} new={cell} old={prev_cell}'
            
        self.total_overlap += overlap
        
        irun, ioff = coord
        run = self.run_at(irun)
        log.debug(f'Assigning coord={coord} overlap={overlap} nc={"X" if prev_cell is None else prev_cell.nchan}->{cell.nchan} old={baselinestr(prev_cell)} requested={baselinestr(cell)}')

        assert run is not None
        assert self.cell_at(coord) is None, f'Expected empty cell at {coord} but got {self.cell_at(coord)}'

        inserted_cell, new_cell = run.insert_cell(ioff, cell)

        return inserted_cell, new_cell
    
    def chan_start_for_new_cell(self, cell):
        '''
        Returns the chan start we should use for a new cell
        Hacky currently
        '''
        cstart = chan_with_best_overlap(cell, self.chan_positions, self.chan_width)
        return cstart
    
    def add_cell(self, cell):
        '''
        Find a run to put the given cell in.
        If none exists, add a new run
        Returns the a tuple (newcell, leftover_cell) - the leftover_cell is always none.
        '''
        # let's just assume prev_plan is defined for now
        assert self.prev_plan is not None
        target_run = None

        for irun, run in enumerate(self.__runs):
            if run is None:
                continue
                
            run_overlap = channel_overlap(cell.chan_start, cell.chan_end, run.chan_start, run.chan_end)
            prev_run = self.prev_plan.run_at(irun)
            
            # TODO - we might consider adding a cell to an existing run if it has less overlap
            # but it probably isn't important
            # for now we require full overlap
            # If this cell overlaps with the run, then let's see if we cean put it in.
            if run_overlap == cell.nchan:            # if we overlap nicely

                target_icell = -1
                if prev_run is None:
                    # this is a new run that we've added to contain new stuff. I think we just keep using it
                    if run.ncell < self.nuvwide:# We can only add to this new cell if it isn't full
                        target_icell = run.next_empty_index
                else:
                    # loop through cells of current and previous run and find the index where both are empty
                    for icell, (c1, c2) in enumerate(zip(prev_run.cells, run.cells)):
                        if c1 is None and c2 is None:
                            target_icell = icell
                            break
                        
                # if there is such an index, then we're good
                # otherwise, we drop out of the IF and keep going
                if target_icell >= 0:
                    target_irun = irun
                    target_run = run 
                    log.debug(f'Found spare cell in irun={target_irun} {prev_run} at icell={target_icell}')
                    break

        # No available run. make a new one.
        if target_run is None:
            cstart = self.chan_start_for_new_cell(cell)
            target_run, target_irun = self.add_run(cstart)
            target_icell = 0
            log.debug('No space for cell. Created new run idx=%d for %s at cstart=%s', target_irun, cell, cstart)
        
        target_coord =  (target_irun, target_icell)
        
        self.set_cell_coord(cell, target_coord)
        
        return cell, None
                
    
    def delete_empty_runs(self):
        '''
        Remove runs that are empty after all the planning
        '''
        for irun, run in enumerate(self.__runs):
            if run is not None and run.ncell == 0:
                self.__runs[irun] = None
            
        
    def get_cell_coord(self, cell:BaselineCell):
        '''
        Returns the coordinate of the given BaselineCell
        Raises ValueError if the BaselineCell is not found
        '''
        coord = None
        for irun, run in enumerate(self.__runs):
            if run is None:
                continue
            if cell in run.cells:
                cellidx = run.cells.index(cell)
                coord = (irun, cellidx)
                
        if coord is None:
            raise ValueError(f'Cell not found {cell}')
            
        assert self.cell_at(coord) == cell
            
        return coord
    
    def delete_cell_coord(self, coord):
        '''
        Doesn't do anything but checks that the coordinate was valid in teh previus run and no longer in the future 
        one
        '''
        irun, ioff = coord
        assert self.cell_at(coord) is None
        if self.prev_plan is not None:
            assert self.prev_plan.cell_at(coord) is not None

    def __str__(self):
        pcoverlap = self.total_overlap / self.nchan * 100
        s = f'FdmtPlanContainer nruns={self.nruns} nc={self.nchan} overlap={self.total_overlap}={pcoverlap:0.1f}%'
        return s

    __repr__ = __str__
            
            
    
def print_cell_diffs(p1cells, p2cells):
    '''
    Debug print function
    '''
    print('BLID', p1cells[0].blid)
    for i in range(max(len(p1cells), len(p2cells))):       
        print(i,rangestr(i, p1cells), rangestr(i, p2cells))
        
    

def migrate_plan(plan1, plan2, blids=None):
    '''
    Migrate from plan1 to plan2
    Assumes plan2 is empty to start
    Tries to match UVWs in plan2 up with plan1
    If there are extra cells it will add them to plan2
    Leftover UVs will be deleted
    '''
    bl1 = plan1.blids
    bl2 = plan2.blids
    if bl1 != bl2:
        log.warning(f'Baselines changed between plans: {len(bl1)} =? {len(bl2)} - d1={bl1 - bl2} d2={bl2 - bl1} d12={bl1 ^ bl2}. Probably due to U=0 V=0 bug')
        
    if blids is None:
        blids = plan1.blids.union(plan2.blids) # set union  - some blids can be in ond and not the other

    for blid in blids:
        p1cells = plan1.uvcells_of(blid)
        p2cells = find_cells(plan2.uvcells, blid)
        log.debug(f'Plan2 nchan={plan2.nchan} overlap={plan2.total_overlap} next run={plan2.next_empty_run_index}')
        #print_cell_diffs(p1cells, p2cells)
        extra_p2_cells = []
        
        # assign as many of plan2 cells to their coordinates as we can in plan1
        while len(p2cells) > 0:
            if len(p1cells) == 0: # nothing left to assign. We'll have leftovers
                break
                
            p2c = max(p2cells, key=lambda c:c.nchan)
            p2cells.remove(p2c)

            p1c = get_cell_with_highest_overlap(p2c, p1cells)            
            overlap = cell_channel_overlap(p1c,p2c)
            #print(f"best overlap of {p2c} with n={len(p1cells)}={p1cells} is {p1c} overlap={overlap}")
            
            if overlap == 0:
                # Remove. there's nothing to overlap. We should save it to add to new runs
                # we need to remove from p2cells otherwise the loop never ends, but add it later
                # so it gets added to new runs.
                extra_p2_cells.append(p2c)
            else:
                p1coord = plan1.get_cell_coord(p1c)
                inserted_p2c, new_p2c = plan2.set_cell_coord(p2c, p1coord)
                p1cells.remove(p1c)
                if new_p2c is not None:
                    p2cells.append(new_p2c)                            
            
        # we might end up with leftover p2cells. In which case we need to add them as new cells to the plan
        p2cells.extend(extra_p2_cells)
        log.debug(f'There are {len(p2cells)} cells that need to be added: {p2cells}')

        while len(p2cells) > 0:            
            p2c = max(p2cells, key=lambda c:c.nchan)
            inserted_p2c, new_p2c = plan2.add_cell(p2c)
            p2cells.remove(p2c)
            if new_p2c is not None:
                p2cells.append(new_p2c)
            
        # for p1 cells that don't have assignments we need to delete them
        log.debug(f'There are {len(p1cells)} cells that need to be deleted: {p1cells}')

        for p1c in p1cells:
            p1coord = plan1.get_cell_coord(p1c)
            plan2.delete_cell_coord(p1coord)
            
            
    plan2.delete_empty_runs()
    log.info('Migrated plan from %s to %s', plan1, plan2)

    if plan1.blids == plan2.blids: # all bets are off if the blids changed
        expected_nchan = sum(c.nchan for c in plan2.uvcells)
        assert plan2.nchan == expected_nchan, f'Didnt plan all cells. Expected={expected_nchan} actual={plan2.nchan}'
            
    return plan2

def migrate_pipeline_plans(plan1, plan2, new_runs):
    log.info('MIgrate plan')
    container1 = FdmtPlanContainer(plan1.fdmt_plan, plan1.fdmt_plan.runs, None)
    container2 = FdmtPlanContainer(plan2.fdmt_plan, new_runs, container1)
    container2 = migrate_plan(container1, container2)
    migrated_runs = container2.valid_runs # only up the last valid index
    return migrated_runs

def migrate_fdmt_plans(plan1, plan2, new_runs):
    log.info('Migrate plan from %s to %s with %d new runs', plan1, plan2, len(new_runs))
    container1 = FdmtPlanContainer(plan1, plan1.runs, None)
    container2 = FdmtPlanContainer(plan2, new_runs, container1)
    container2 = migrate_plan(container1, container2)
    #migrated_runs = container2.valid_runs
    # shorten to make it only as long as it needs to be
    return container1, container2


class FdmtRun:
    '''
    FDMT run - Maintains a list of BaselineCell that will dedipsersed together. cells NUVWIDE
    They all start blank.
    You can append to them with add_cell
    Or set with set_cell(idx, cell)
    If it's already full, or that cell is set then its an exception
    '''
    def __init__(self, plan, chan_start):
        self.plan = plan
        self.nuvwide = self.plan.nuvwide
        self.nf = self.plan.pipeline_plan.nf
        self.nuvmax = self.plan.nuvmax

        self.__cells = [None for i in range(self.plan.nuvwide)]
        assert chan_start >= 0 and chan_start < self.plan.pipeline_plan.nf, f'Invalid chan_start {chan_start}'
        self.chan_start = chan_start
        assert self.max_idm <= plan.pipeline_plan.ndout, f'NDOUT ={plan.pipeline_plan.ndout} is too small - needs to be at least {self.max_idm} for fch1={self.fch1} fmin={plan.pipeline_plan.fmin} fmax={plan.pipeline_plan.fmax} dmax={self.plan.pipeline_plan.dmax}'


    @classmethod
    def from_cells(self, cells, plan, chan_start=None):
        if chan_start is None:
            mincell = min(cells, key=lambda cell:cell.chan_start)
            chan_start = mincell.chan_start
            
        run = FdmtRun(plan, chan_start)
        for c in cells:
            assert len(cells) <= plan.nuvwide, f'Too many cells: nuvwide={plan.nuvwide} ncells={len(cells)}'
            run.add_cell(c)

        return run

    @property
    def cells(self):
        '''
        Returns a copy of the cells array so noone fiddles with it
        Length will be NUVWIDE and may contain None
        '''
        return list(self.__cells)

    @property
    def defined_cells(self):
        '''
        Returns an iterator of the cells without the None values
        '''
        return filter(lambda c: c is not None, self.__cells)

    def append_cell(self, cell):
        '''
        Adds the cell after the last valid cell (probably useless)
        '''
        idx = -1
        for i, c in enumerate(self.__cells):
            if c is not None:
                idx = i

        if idx == len(self.__cells) - 1:
            raise ValueError(f"FdmtRun is already full! idx={idx} len={len(self.__cells)} {self.__cells}")
   
        self.set_cell(idx+1, cell)
        assert self.ncell <= self.plan.nuvwide
        return self

    @property
    def next_empty_index(self):
        cellidx = -1
        for icell, cell in enumerate(self.__cells):
            if cell is None:
                cellidx = icell
                break

        assert cellidx != -1
        return cellidx

    def add_cell(self, cell):
        '''
        Adds the given cell to the first available slot

        Returns the cell index where it was inserted
        '''
        assert self.ncell < self.nuvwide, 'This run is full!'
        cellidx = self.next_empty_index
        self.set_cell(cellidx, cell)
        return cellidx

    def set_cell(self, idx, cell):
        '''
        Set a baseline cell at the given index
        rasies exception if something's already there
        Raises exception if the cell channel range exceeds the run channel range
        '''

        assert cell.chan_start >= self.chan_start
        assert cell.chan_end <= self.chan_end, f'End channel out of range: {cell} this={self.chan_start}-{self.chan_end}'
        assert cell.nchan <= self.plan.pipeline_plan.ncin, f'Input cell has too many channels: {cell} ncin={self.plan.pipeline_plan.ncin}'

        if self.__cells[idx] is not None:
            raise ValueError(f'FdmtRun Cell at {idx} is already set to {self.__cells[idx]} trying to set with {cell}')
        
        self.__cells[idx] = cell
        cell.fdmt_run = self

        return self


    def insert_cell(self, idx, cell):
        '''
        Sets a baseline cell at the given index.
        The cell must have a chan_start >= the chan_start of this FDMT run
        If the cell end channel exceeds the chan_end of this run, the input cell is split and set at the given index
        and the remaining cell is returned.
        If the cell is not split then None is returned
        '''
        assert cell.chan_start >= self.chan_start, f'Invalid chan start {self} {cell}'
        new_cell = None
        if cell.chan_end > self.chan_end:
            insert_cell, new_cell = cell.split_at(self.chan_end)
        else:
            insert_cell = cell

        self.set_cell(idx, insert_cell)
        
        return insert_cell, new_cell

    @property
    def chan_end(self):
        '''
        Final channel, inclusive
        '''
        return self.chan_start + self.plan.pipeline_plan.ncin - 1

    @property
    def fch1(self):
        return self.plan.pipeline_plan.freqs[self.chan_start]

    @property
    def total_overlap(self):
        over = 0
        for uv in self.defined_cells:
            overlap = calc_overlap(uv, self.chan_start,self.plan.pipeline_plan.ncin)
            #log.debug('Cell chan_start %s %s %s-%s overlap=%d', self.chan_start, uv, uv.chan_start, uv.chan_end, overlap)
            assert overlap > 0
            over += overlap
        return over
        

    @property
    def ncell(self):
        # count the number of nonzero values in cells
        return sum(c is not None for c in self.cells)
    
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
        ncells = self.ncell
        return 'ncells={ncells} fch1={self.fch1} chan_start={self.chan_start} total_overlap={self.total_overlap}'.format(self=self, ncells=ncells)

    __repr__ = __str__

class FdmtPlan(object):
    def __init__(self, uvcells, pipeline_plan, prev_pipeline_plan=None):
        self.pipeline_plan = pipeline_plan
        nuvwide = self.pipeline_plan.nuvwide
        ncin = self.pipeline_plan.ncin
        self.nuvwide = nuvwide
        self.ndout = self.pipeline_plan.ndout
        self.ncin = ncin
        self.nuvmax = self.pipeline_plan.nuvmax
        self.uvcells = uvcells
        self.prev_pipeline_plan = prev_pipeline_plan
        self.zero_cell = None
        self.runs = []
        
        #runs = greedy_create_fdmt_runs(self, uvcells, nuvwide, ncin)
        runs = fixed_create_fdmt_runs(pipeline_plan, self, uvcells, nuvwide, ncin) # even for migrated data we need initial runs

        log.info('Creating FdmtPlan for %d cells plan=%s prev_plan=%s We have %d initial runs', len(uvcells), pipeline_plan, prev_pipeline_plan, len(runs))
        
        if prev_pipeline_plan is None:
            self.set_runs(runs)
        else:
            self.container1, self.container2 = migrate_fdmt_plans(prev_pipeline_plan.fdmt_plan, self, runs)
            migrated_runs = self.container2.valid_runs
            self.set_runs(migrated_runs)

    def __str__(self):
        return f'FdmtPlan ncell={len(self.uvcells)} zero cell={self.zero_cell} nruns={len(self.runs)}'

    __repr__ = __str__
        

    def set_runs(self, runs):
        uvcells = self.uvcells
        nuvwide = self.nuvwide
        ncin = self.ncin
        nuvmax = self.nuvmax
        nuvrest_max = self.pipeline_plan.nuvrest_max

        self.nuvrest_valid = len(runs)
        self.runs = runs
        
        # find a cell with zero in it
        self.zero_cell = None
        for irun, run in enumerate(runs):
            if run is None:
                continue
            
            if run.ncell < nuvwide - 1:
                self.zero_cell = (irun, run.next_empty_index)
                log.info('Found zero cell among runs: %s', self.zero_cell)
                break

        # If we still haven't foudn another UVCell, we need to add another FDMT run
        # Which seems an aweful waste, but go figure.
        if self.zero_cell is None: 
            runs.append(FdmtRun(self, 0)) # Add FDmt Run
            self.zero_cell = (len(runs)-1, 0)
            log.info('Had to add new run just for zero cell: %s', self.zero_cell)
            
        log.info("FDMT zero cell is %s=%s", self.zero_cell, self.zero_cell[0]*nuvwide+self.zero_cell[1])

        assert self.zero_cell != None
        assert self.zero_cell[0] < self.pipeline_plan.nuvrest_max, 'Not enough room for FDMT zero cell'
        assert self.zero_cell[1] < nuvwide
        assert self.zero_cell[0] < self.nruns, f'Zero Cell unexpected {self.zero_cell}[0] should be less than {self.nruns}'
        assert self.zero_cell[1] < ncin
        assert runs[self.zero_cell[0]].cells[self.zero_cell[1]] is None
        assert self.get_cell(self.zero_cell) is None

        nruns = self.nruns
        nuvtotal = nruns*nuvwide
        ndout = self.pipeline_plan.ndout
        nd = self.pipeline_plan.nd
        nt = self.pipeline_plan.nt
        #square_history_size = ndout*nuvtotal*(nt + nd)
        valid_runs = [run for run in runs if run is not None]
        square_history_size = sum(nuvwide*(nd + nt)*ndout for run in runs)
        minimal_history_size = sum(nuvwide*(run.max_offset+ nt)*run.max_idm for run in valid_runs)
        efficiency = float(len(uvcells))/float(nuvtotal)
        required_efficiency = float(nuvtotal)/8192.0
        
        log.info('FDMT plan has ntotal=%d of %d runs with packing efficiency %f. Grid read requires efficiency of > %f of NUV=8192. History size square=%d minimal=%d =%d 256MB HBM banks', nuvtotal, nruns, efficiency, required_efficiency, square_history_size, minimal_history_size, minimal_history_size*4/256/1024/1024)
        freqs = [self.pipeline_plan.fmin if run is None else run.fch1 for run in runs]
        self.fdmt_lut = self.make_fdmt_lut(freqs)
        self.nuvtotal = nuvtotal
        self.total_nuvcells = sum([run.ncell for run in valid_runs])

        # create UVMAP for more rapid lookups 
        uvmap = {}
        for irun, run in enumerate(self.runs):
            if run is None:
                continue
            
            for icell, cell in enumerate(run.cells):
                if cell is None:
                    continue
                irc = (irun, icell)
                uvmap[cell.uvpix_upper] = uvmap.get(cell.uvpix_upper, [])
                uvmap[cell.uvpix_upper].append(irc)
                assert irc is not None

        self.__uvmap = uvmap

    @property
    def nruns(self):
        return len(self.runs)

    @property
    def default_fdmt_lut(self):
        '''
        Returns an FDMT LUT that can be programmed into the hardware for the at the
        minimum frequency of the plan. Setup for NUVREST_MAX = all runs
        Useful for clearing memory on the card
        '''
        freqs = np.ones(self.pipeline_plan.nuvrest_max)*self.pipeline_plan.fmin
        lut = self.make_fdmt_lut(freqs)
        return lut


    def make_fdmt_lut(self, freqs):
        '''
        Make and FDMT LUT for the given frequencies
        '''
        lut = make_fdmt_lut(freqs, self.pipeline_plan.foff, self.ncin, self.ndout, self.nuvwide)
        return lut


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
        cell_coords2 = self.__uvmap.get(uvpix)
        
        #cell_coords = []
        #for irun, run in enumerate(self.runs):
        #    for icell, cell in enumerate(run.cells):
        #        if cell.uvpix_upper == uvpix:
        #            cell_coords.append((irun, icell))


        #print(uvpix, 'Version1', cell_coords, 'Verion2', cell_coords2)
        #assert cell_coords == cell_coords2
        #assert cell_coords2 is not None, f'No coordinates in UVMAP at uvpix={uvpix}'
        if cell_coords2 is None:
            log.warning('Could not find UV for uvpix: %s', uvpix)

        return cell_coords2

    def get_cell(self, cell_coord):
        '''
        Returns baseline cell of given FDMT run and cell cell coord (run, cell)
        if cell coord is the zero cell, it returns None
        '''
        irun, icell = cell_coord
        if cell_coord == self.zero_cell:
            return None
        else:
            return self.runs[irun].cells[icell]




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
