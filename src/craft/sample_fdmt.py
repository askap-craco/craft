#!/usr/bin/env python
"""
Classes that implement the sample FDMT

Copyright (C) CSIRO 2020
"""
from abc import ABCMeta, abstractmethod, abstractproperty
import numpy as np
import logging
import math

__author__ = "Keith Bannister <keith.bannister@csiro.au>"
def is_power_of_2(x):
    return math.log(x, 2).is_integer()

    
class SampleFdmt(object, metaclass=ABCMeta):
    '''
    Abstract base class (ABC) representing  a sample based FDMT based on the supplied block-based FDMT.

    The sample FDMT internally has the concept of FIFOs in it. Subclasses must implement their own scheme for
    creating a buffer for storing the FIFO and shifting() and reading() from the FIFOS. 

    This class implements 2 methods:
    - fdmt_process() - processes a single time sample and produces the dedispersed output
    - execute() - which takes in a block of data, initialises it, and runs fdmt_process on each time sample

    '''
    
    def __init__(self, thefdmt):
        '''
        Creates the class
        
        :thefdmt: block based FMDT  which contains much of the logic 
        to compute all the offsets and has all the numbers in it
        '''
        self.thefdmt = thefdmt
        if not is_power_of_2(self.thefdmt.n_f):
            raise NotImplementedError("Don't support non-power-of-2 number of channels yet")
        
        
    @abstractmethod
    def shift(self, iterno, d, c, value):
        '''
        Take the supplied value and push it into a FIFO indexed by interno, d, and c at t=0. The old value drops out.
        
        :iterno: Iteration number (int)
        :d: DM value (int)
        :c: channel (int)
        :value: value to push on at t = 0
        '''
        pass
    
    @abstractmethod
    def read(self, iterno, d, c, t):
        '''
        Read the value of the FIFO
        
        :iterno: Iteration value to read from (int)
        :d: DM value (int)
        :c: Channel (int)
        :t: Time offset to read (int). t=0 is the most recently shifted value
        '''
        pass

    @abstractproperty
    def buffer_size(self):
        '''
        Returns the number of entries in all the buffers allocated to this SampleFdmt
        '''
        pass
    
    def fdmt_process(self, din):
        '''
        Takes in a single initialised time sample and returns the FDMT of it.
        
        :din: np.array of shape (NCHAN, ND_IN)
        :dout: np array of shape (ND)
        '''
        thefdmt = self.thefdmt

        assert din.shape[0] == thefdmt.n_f
        assert din.shape[1] == thefdmt.init_delta_t
        
        # Push the input data in to the iteration 0 FIFOs
        for c in range(din.shape[0]):
            for d in range(din.shape[1]):
                self.shift(0, d, c, din[c, d])

        niter = len(thefdmt.hist_nf_data)
        dout = np.zeros(thefdmt.max_dt)

        for iterno, theiter in enumerate(thefdmt.hist_nf_data):
            for output_channel in range(len(theiter)):
                chanconfig = thefdmt.hist_nf_data[iterno][output_channel][-1]
                for out_d, config in enumerate(chanconfig):
                    # id1, id2, offset are values in the lookup table
                    in_d1 = config[1]
                    in_d2 = config[3]
                    time_offset = config[2]
                    in_chan1 = 2*output_channel
                    in_chan2 = 2*output_channel+1
                    
                    # Read values from FIFOs at d, c read from 
                    v1 = self.read(iterno, in_d1, in_chan1, 0)
                    v2 = self.read(iterno, in_d2, in_chan2, time_offset)
                    vout = v1 + v2
                    if iterno == niter - 1: # final iteration write to output
                        dout[out_d] = vout
                    else:
                        self.shift(iterno+1, out_d, output_channel, vout) 

        return dout


    def execute(self, inblock):
        '''
        Takes in a block of data, initialises it and produces a block of data by running fdmt_process on it for each
        time sample

        :inblock: Block of data np.array of shape (NCHAN, NT)
        :returns: Block of data - np.array of shape (ND, NT)
        '''
        
        thefdmt = self.thefdmt
        dinit = thefdmt.initialise(inblock)
        # Dinit shape is nf, nd, nt
        ntout = dinit.shape[2]
        out = np.zeros((thefdmt.max_dt, ntout))
        for t in range(ntout):
            out[:, t] = self.fdmt_process(dinit[:, :, t])

        return out
    
    __call__ = execute

    
class MaxFifo(SampleFdmt):
    '''
    Makes 1 giant buffer [NITER, ND, NC, MAX_FIFO_SIZE]. 
    Evalent to a FIFO for every iteration, channel and DM equal to the largest FIFO.
    Its *HUGE*
    '''
    def __init__(self, thefdmt):
        super(MaxFifo, self).__init__(thefdmt)
        nd = thefdmt.max_dt
        nc = thefdmt.n_f
        max_offset = 0
        niter = len(thefdmt.hist_nf_data)

        for curr_iterno, theiter in enumerate(thefdmt.hist_nf_data):
            ochan = 0 # Largest offset will be for the 0 channel
            chanconfig = thefdmt.hist_nf_data[curr_iterno][ochan][-1]
            offsets = np.array([config[2] for config in chanconfig]) # offset vs IDT
            max_offset = max(max_offset, offsets.max())

        # Gosh that's a big buffer
        max_fifo_size = max_offset + 1
        self.buffer = np.zeros((niter, nd, nc, max_fifo_size))
        
    def shift(self, iterno, d, c, value):
        self.buffer[iterno,d,c,1:] = self.buffer[iterno,d,c,:-1]
        self.buffer[iterno,d,c,0] = value
    
    def read(self, iterno, d, c, t):
        # Look up value in a giant array
        value = self.buffer[iterno, d, c, t]
        return value
    
    def buffer_size(self):
        return self.buffer.size
    

class MaxFifoPerIteration(SampleFdmt):
    '''
    Makes 1 buffer per iteration. Each FIFO for a given iteration is the size of the largest FIFO in that iteration
    '''
    def __init__(self, thefdmt):        
        super(MaxFifoPerIteration, self).__init__(thefdmt)        
        self.buffers = []        
        self.fifo_size = 0
        for curr_iterno, theiter in enumerate(thefdmt.hist_nf_data):
            ochan = 0
            chanconfig = thefdmt.hist_nf_data[curr_iterno][ochan][-1]
            offsets = np.array([config[2] for config in chanconfig]) # offset vs IDT
            fifo_sizes = offsets + 1
            self.fifo_size += fifo_sizes.sum()
            # TODO: Find some way of predicting the offset vs iteration give fdmt input paramters
            # equal for all channels
            state_shape = thefdmt.hist_state_shape[curr_iterno]
            nc = state_shape[0] # number of output channels in this iteration
            nd = state_shape[1] # maximum number of DMs in this iteration
            buf = self.make_buffer(nd, nc, fifo_sizes)            
            print('Iteration', curr_iterno, 'buffer ', buf.shape, 'size', buf.size, ' maxbuf', fifo_sizes.max(), 'state_shape', state_shape)
            self.buffers.append(buf)
            
    def buffer_size(self):
        return sum([np.array(b.shape).prod() for b in self.buffers])
    
    def make_buffer(self, nd, nc, fifo_sizes):
        max_fifo_size = fifo_sizes.max()
        assert max_fifo_size == fifo_sizes[-1]
        buf = np.zeros((nd, nc, max_fifo_size+1)) # OK - let's make it a square for now. 
        return buf
        
    def shift(self, iterno, d, c, value):
        # Shift everything by along the time axis
        self.buffers[iterno][d, c, 1:] = self.buffers[iterno][d, c, 0:-1]
        self.buffers[iterno][d, c, 0] = value
        
    def read(self, iterno, d, c, t):
        v = self.buffers[iterno][d,c,t]
        return v
    
class IndividualFifos(SampleFdmt):
    '''
    Makes 1 buffer per FIFO of exactly the right length for the output of each node
    This will only make 1 output FIFO per node that is the correct size
    Usually a FIFO will have a read at t=0 and a read at t=FIFO_LENGTH-1, but 
    occasionally it will have multiple reads at 0 < t < FIFO_LENGTH
    '''
    def __init__(self, thefdmt):
        super(IndividualFifos, self).__init__(thefdmt)        
        self.fifos = {} # Dictionary of FIFOS (cheater) key=(iterno, d, c) 
        self.__buffer_size = 0
        for curr_iterno, theiter in enumerate(thefdmt.hist_nf_data):
            for output_channel in range(len(theiter)):
                chanconfig = thefdmt.hist_nf_data[curr_iterno][output_channel][-1]
                for idt, config in enumerate(chanconfig):
                    in_d1 = config[1]
                    in_d2 = config[3]
                    time_offset = config[2]
                    fifo_size = time_offset + 1
                    in_chan1 = 2*output_channel
                    in_chan2 = 2*output_channel+1
                    # you only ever read t=0 sample from inchan1
                    self._create_fifo(curr_iterno, in_d1, in_chan1, 1)
                    self._create_fifo(curr_iterno, in_d2, in_chan2, fifo_size)

    def _create_fifo(self, iterno, d, c, size):
        '''
        Makes a fifo of the requested size and puts it in the fifo dictionary
        and keeps track of the total buffer size
        If the fifo already exists in the dictionary, it is resized to provide enough room for the maximum size
        '''
        key = (iterno, d, c)
        if key in self.fifos:
            fifo = self.fifos[key]
            if len(fifo) < size:
                # Remove the old FIFO and we'll make a new one
                self.__buffer_size -= len(fifo)
                del self.fifos[key]
                fifo = np.zeros(size)
                self.fifos[key] = fifo
                self.__buffer_size += len(fifo)
            else: # existing FIFO OK size
                pass
        else: # no FIFO for that key yet
            fifo = np.zeros(size)
            self.fifos[key] = fifo
            self.__buffer_size += len(fifo)
    
    def _get_fifo(self, iterno, d, c):
        # Need to do a 3 pointer derefrences to find the FIFO of interest
        # oR just lookup in a dictionary for now
        fifo = self.fifos[(iterno, d, c)]
        return fifo

    def buffer_size(self):
        return self.__buffer_size
        
    def shift(self, iterno, d, c, value):
        try:
            fifo = self._get_fifo(iterno, d, c)
            fifo[1:] = fifo[:-1]
            fifo[0] = value
        except KeyError: # Happens occasionally. Not all iterno/d/c have a fifo
            pass
            
    def read(self, iterno, d, c, t):
        try:
            fifo = self._get_fifo(iterno, d, c)
            v = fifo[t]
        except:
            raise ValueError('Could not read value for iterno={} d={} c={} t={} '.format(iterno, d, c, t))
        return v

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    

if __name__ == '__main__':
    _main()
