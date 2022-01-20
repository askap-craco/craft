#!/usr/bin/env python
"""
Generates HLS code for an FDMT. based on 'FDMT one sample at a time.ipynb'

A work in progress. It's a mess but will hopefully improve with time

Copyright (C) CSIRO 2020
"""
from pylab import *
import matplotlib as mpl
from collections import OrderedDict
import logging
import os
import numpy as np
from scipy import constants
from . import fdmt # you'll need to have ../python in  you PYTHONPATH
from graphviz import Digraph
import datetime

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

def makedirs(d, exist_ok=False):
    if os.path.exists(d) and not exist_ok:
        raise ValueError('Path already exists')

    try:
        os.makedirs(d)
    except: # ignore error anyway
        pass

    

def mkline(l, indent, endl='\n'):
    s = ''
    if isinstance(l, Cblock):
        s = str(l)
    else:
        s = ' '*indent + str(l)
        if not s.endswith(endl):
            s += endl

    return s

class Cblock():
    def __init__(self, parent, hdr, indent=0, formatdict=None):
        if formatdict:
            self.hdr = hdr.format(**formatdict)
        else:
            self.hdr = hdr

        self.lines = []
        self.parent = parent
        self.indent = indent


    def __str__(self):
        s = mkline(self.hdr, self.indent)
        s += mkline('{', self.indent)
        for l in self.lines:
            s += mkline(l, self.indent + 4)
        s += mkline('}', self.indent)
        return s

    def __iadd__(self, line):
        self.lines.append(line);

        return self

    def block(self, hdr):
        blk = Cblock(self, hdr, indent=self.indent + 4)
        self.lines.append(blk)
        return blk

class Code:
    def __init__(self, preamble='', hfile='FDMT_PROCESS'):
        self.functions = []
        self.funcdefs = []
        self.cfile = ''
        self.preamble = preamble
        self.hfile = hfile

    @property
    def header(self):
        body = ''.join(mkline(func.hdr, indent=0, endl=';\n') for func in self.functions)
        if self.hfile:
            s = '''
#ifndef _{hfile}
#define _{hfile}

{body}
#endif
'''.format(hfile=self.hfile, body=body)
        else:
            s = body

        return s

    @property
    def code(self):
        s = self.preamble + '\n\n'.join(map(str, self.functions))
        return s

    def cfunction(self, funcdef, formatdict=None):
        thefunc = Cblock(self, funcdef, indent=0, formatdict=formatdict)
        self.functions.append(thefunc)
        return thefunc

def fmt(i,c,d):
    return 'I{}D{}C{}'.format(i,d,c)

ignore_ant = ['ak31','ak32','ak33','ak34','ak35','ak36']

# List of configuraitons to dump by Default
# (Nd, Nchan, )
configs = ((2,2), (4, 4), (8,8), (16,16), (32,32), (64,64),(128,128),(256,256),(1024,256),(1024,288),(256,288),(32,288), (512, 256), (1024, 32), (1024, 64), (1024, 128))

class FdmtDag(object):
    '''
    A directed acyclic graph version of the suplied fdmt
    '''

    def __init__(self, thefdmt):
        self.thefdmt = thefdmt
        maindot = Digraph(comment='FDMT', format='png')
        self.maindot = maindot
        all_offsets = []
        self.all_offsets = all_offsets
        nnodes = 0
        nplotnodes = 0
        output_fifo_sizes = {}
        self.output_fifo_sizes = output_fifo_sizes

        # make initial nodes
        ishape = thefdmt.hist_state_shape[0]
        dot = Digraph(name='cluster_iter0')
        dot.attr(label='Iteration 0')
        maindot.subgraph(dot)
        for c in range(ishape[0]):
            for d in range(ishape[1]):
                maindot.node(fmt(0, c, d))
                nnodes += 1
                nplotnodes += 1


        for iterno, nfd in enumerate(thefdmt.hist_nf_data):
            out_shape = thefdmt.hist_state_shape[iterno+1]
            in_shape = thefdmt.hist_state_shape[iterno]
            nchan, ndt, nt_out = out_shape
            #print 'Iteration {} in={}={} out={}={}'.format(iterno, in_shape, in_shape[0:2].prod(), out_shape, out_shape[:2].prod())
            dot = Digraph(name='cluster_iter{}'.format(iterno+1))
            dot.attr(label='Iteration {}'.format(iterno+1))
            maindot.subgraph(dot)

            for ochan in range(nchan):
                chanconfig = thefdmt.hist_nf_data[iterno][ochan][-1]
                #print '\tOut channel {}'.format(ochan)
                last_id1 = -1
                last_id2 = -1
                for idt, config in enumerate(chanconfig):
                    _, id1, offset, id2, _, _, _ = config
                    do_copy = id2 == -1
                    inchan1 = 2*ochan
                    inchan2 = inchan1+1
                    id1_hit = ''
                    id2_hit = ''
                    if last_id1 == id1:
                        id1_hit = '*'
                    if last_id2 == id2:
                        id2_hist = '*'

                    all_offsets.append(offset)
                    style = 'dashed' if do_copy else 'solid'
                    nnodes += 1
                    label = None if offset == 0 else 'd={}'.format(offset)
                    n1 = fmt(iterno, inchan1, id1)
                    n2 = fmt(iterno, inchan2, id2)
                    nout = fmt(iterno+1,ochan, idt)

                    # For a given sum - only maintain 1 fifo of the given max length
                    ff = output_fifo_sizes.get(n2, [])
                    ff.append(offset)
                    output_fifo_sizes[n2] = ff

                    maindot.node(fmt(iterno+1,ochan, idt), style=style)
                    nplotnodes += 1
                    maindot.edge(n1, nout)
                    if not do_copy:
                        color = 'black' if offset == 0 else 'red'
                        maindot.edge(n2, nout,label=label, color=color)

                    last_id1 = id1
                    last_id2 = id2


    @property
    def nfifo_outputs(self):
        return {k:max(ff) for k, ff in output_fifo_sizes.items()}

    @property
    def bulk_fifo_sizes(self):
        return {k:min(ff) for k, ff in output_fifo_sizes.items()}

    @property
    def fanout_fifo_sizes(self):
        return {k:max(ff) - min(ff) for k, ff in output_fifo_sizes.items()}

    @property
    def num_outputs(self):
        return [len(ff) for ff in list(self.output_fifo_sizes.values())]

    @property
    def max_ff_length(self):
        return [max(ff) for ff in list(self.output_fifo_sizes.values())]

    @property
    def ff_length_range(self):
        return [max(ff) - min(ff) for ff in list(self.output_fifo_sizes.values())]

    @property
    def total_offsets(self):
        return np.array(list(self.nfifo_outputs.values()))

    @property
    def bulk_sizes(self):
        return np.array(list(self.bulk_fifo_sizes.values()))

    @property
    def fanout_sizes(self):
        return np.array(list(self.fanout_fifo_sizes.values()))

    def print_stats(self):
        output_fifo_sizes = self.output_fifo_sizes
        nfifo_outputs = self.nfifo_outputs
        bulk_fifo_sizes = self.bulk_fifo_sizes
        fanout_fifo_sizes = self.fanout_fifo_sizes
        all_offsets = self.all_offsets
        print('Total offsets', sum(all_offsets), 'nfifo outputs', sum(nfifo_outputs.values()), 'largest fifo', max(nfifo_outputs.values()), 'total nodes', nnodes)


    def print_sr_stats(self, all_offsets):
        all_offsets = self.all_offsets
        print("num FIFOs", len(all_offsets))
        print('Total SR entries', sum(all_offsets))
        print('Number of SR length ==0', sum(all_offsets == 0))
        print('Number of SR length ==1', sum(all_offsets == 1))
        print('Number of SR length ==2', sum(all_offsets == 2))
        print('Number of SR length >=1', sum(all_offsets >= 1))
        print('Number of SR length <=16', sum(all_offsets <= 16))
        print('Number of SR length <=32', sum(all_offsets <= 32))
        print('Number of SR length <=64', sum(all_offsets <= 64))
        print('Number of SR length >64', sum(all_offsets >64))
        print('Number of SR length >128', sum(all_offsets >128))
        print('Number of SR length >256', sum(all_offsets >256))
        print('Number of SR length >512', sum(all_offsets >512))
        print('Max SR length', max(all_offsets))


    def print_transfer_stats(self):
        # TODO: Fix this funciton. it wont' work
        transfer_size_bits = 512 # DRAM bus width
        entry_bits = 8 # Bits per entry
        entries_per_transfer = float(transfer_size_bits / entry_bits)
        clocks_per_block = float(nt*max(nchan, nd)/entries_per_transfer) # max(nchan, nd) - is depending on whether it's input bound or output-bound
        transfer_clocks = sum(total_offsets)/entries_per_transfer
        longest_fifo = max(total_offsets)
        worst_case_load_nclk = float(max(longest_fifo, transfer_clocks))
        best_case_load_nclk = float(min(longest_fifo, transfer_clocks))
        nfifos_by_size, fsize= np.histogram(total_offsets, bins=np.arange(0, max(total_offsets) + 1))

        print('Entries per transfer', entries_per_transfer, 'processing clocks_per_block', clocks_per_block)
        print('longest_fifo', longest_fifo, 'transfer_clocks', transfer_clocks, 'worst case load nclks', worst_case_load_nclk)
        print('Worst case Processing efficiency', clocks_per_block/(worst_case_load_nclk + clocks_per_block))
        print('Best possible efficiency', clocks_per_block/(best_case_load_nclk + clocks_per_block))
        print('Entry cache Num BRAMS', sum(total_offsets)*8/18e3*2) # 1 for input + 1 for output
        print('Number of entries that can be loaded from memory while clocking largest fifo', longest_fifo*entries_per_transfer, '=', longest_fifo*entries_per_transfer/float(sum(total_offsets))*100, '%')

    def print_all_sr_stats(self):
        print('*'*8, 'total')
        self.print_sr_stats(self.total_offsets)

        print('*'*8, 'bulk')
        self.print_sr_stats(self.bulk_sizes)

        print('*'*8, 'fanout')
        self.print_sr_stats(self.fanout_sizes)

        print('Total offsets', sum(self.total_offsets), ' total sum operations', len(self.all_offsets))


    def plot_fifo_sizes(self):
        x = hist(self.total_offsets, np.arange(0, max(self.total_offsets), 1), log=True, histtype='step')
        xlabel('Total shift register size')
        ylabel('Number of SRs')

        hist(self.total_offsets, np.arange(0, max(self.total_offsets), 1), log=False, cumulative=True, histtype='step')
        xlabel('FIFO size (elements)')
        ylabel('Cumulative number of FIFOs')


        hist(self.total_offsets, np.arange(0, max(self.total_offsets), 1), log=True, cumulative=False, histtype='step')
        xlabel('FIFO size (elements)')
        ylabel('Number of FIFOs')

        plot(x[1][0:-1], x[0]*(x[1][0:-1]))
        xlabel('FIFO size (elements)')
        ylabel('Total number of elements in fifos of this size')

        plot(x[1][0:-1], np.cumsum(x[0]*(x[1][0:-1])))
        xlabel('FIFO size (elements)')
        ylabel('Total number of elements in fifos < this size')

        hist(num_outputs, np.arange(0, max(num_outputs)+3) - 0.5, log=True)
        xlabel('Number of FIFO outputs')

        hist(max_ff_length, np.arange(0, max(max_ff_length)+3) - 0.5, log=True)
        xlabel('Maximum fifo length')


        hist(ff_length_range, np.arange(0, max(ff_length_range)+3) - 0.5, log=True)
        xlabel('Number of FIFO taps')


class FdmtDagFileIter3(object):
    '''
    Splits the caches up into 16-BRAMs wide = versions
    Returns (hstring, cstring) - header file and c file contents as strings

    Based strongly on the original ipython notebook but it's mess and too
    much to fix right now.

    '''

    def __init__(self, fdmt_dag, fifos_per_group=64, max_cache_depth=512, comment_constants=False, compound_threshold=256, bram_depth=1024):
        self.fdmt_dag = fdmt_dag
        self.thefdmt = fdmt_dag.thefdmt
        self.fifos_per_group = fifos_per_group
        self.max_cache_depth = max_cache_depth
        self.comment_constants = comment_constants
        thefdmt = self.thefdmt
        all_offsets = []
        output_fifo_sizes = OrderedDict()
        now = now=datetime.datetime.now().isoformat()
        cname=type(self).__name__
        hfile = '' # Header file
        # make initial nodes
        preamble = '''
// FDMT produced by hls_fdmt.py on {now}
// Class {cname}
// nd={f.max_dt} nf={f.n_f} fmin={f.f_min} nchan={f.n_f} df={f.d_f} bw={f.bw}
// fifos_per_group={fifos_per_group}
// max_cache_depth={max_cache_depth}
// comment_constants={comment_constants}
'''.format(f=thefdmt, **locals())

        ishape = thefdmt.hist_state_shape[0]
        iters = ''
        header_iterdecl = ''
        queuedecl = ''

        for iterno, nfd in enumerate(thefdmt.hist_nf_data):
            out_shape = thefdmt.hist_state_shape[iterno+1]
            in_shape = thefdmt.hist_state_shape[iterno]
            nchan, ndt, nt_out = out_shape
            #print 'Iteration {} in={} size={} out={} size={}'.format(iterno+1, in_shape[0:2], in_shape[0:2].prod(), out_shape[0:2], out_shape[:2].prod())
            ncout, ndout = out_shape[0:2]
            ncin, ndin = in_shape[0:2]
            sums_done = set()
            iterstart = ''
            iterstart += '//Iteration {iterno}\n'.format(**locals())
            if ncout == 1:
                iterdecl = 'void iteration{iterno}(const fdmt_t in[{ndin}][{ncin}], fdmt_t out[{ndout}])'.format(**locals())
            else:
                iterdecl = 'void iteration{iterno}(const fdmt_t in[{ndin}][{ncin}], fdmt_t out[{ndout}][{ncout}])'.format(**locals())

            header_iterdecl += iterdecl + ';\n'
            iterstart += iterdecl + '\n'
            iterstart += '''{
#pragma HLS inline
'''

            queuepush = '// FIFO push statements\n\n'
            do_sums = ''
            read = ' '*4 + '// Read inputs\n'
            for c in range(ncin):
                for d in range(ndin):
                    read += '    fdmt_t {} = in[{}][{}];\n'.format(fmt(iterno, c, d), d,c)

            for ochan in range(nchan):
                chanconfig = thefdmt.hist_nf_data[iterno][ochan][-1]
                #print '\tOut channel {}'.format(ochan)
                last_id1 = -1
                last_id2 = -1

                do_sums += '\n // Output channel {}\n'.format(ochan)
                for idt, config in enumerate(chanconfig):
                    _, id1, offset, id2, _, _, _ = config
                    do_copy = id2 == -1
                    inchan1 = 2*ochan
                    inchan2 = inchan1+1
                    all_offsets.append(offset)

                    n1 = fmt(iterno, inchan1, id1)
                    n2 = fmt(iterno, inchan2, id2)
                    nout = fmt(iterno+1,ochan, idt)

                    # For a given sum - only maintain 1 fifo of the given max length
                    ff = output_fifo_sizes.get(n2, [])
                    ff.append(offset)
                    output_fifo_sizes[n2] = ff

                    if do_copy:
                        do_sums += 'fdmt_t {} = {};\n'.format(nout, n1)
                    else:
                        if offset == 0:
                            do_sums += '    fdmt_t {} = {} + {};\n'.format(nout, n1, n2)
                        else:
                            do_sums += '    fdmt_t {} = {} + {}_fifo.read({});\n'.format(nout, n1, n2, offset-1)

                    sums_done.add((ochan, idt))


            # Find FIFOS for this iteration
            myfifos = [f for f in output_fifo_sizes if f.startswith('I{}'.format(iterno))]
            for infmt in myfifos:
                #ff_sizes = output_fifo_sizes[infmt]
                #queuedecl += 'static fdmt_fifo<{},{}> {}_fifo;\n'.format(min(ff_sizes), max(ff_sizes), infmt)
                #queuedecl += '#pragma HLS ARRAY_PARTITION variable={}_fifo dim=1 complete\n'.format(infmt)
                # add 1 because of how lengths are defined. Ap_shift_reg wants size, fdmt_fifo wanted max_read_value
                queuepush += '    {}_fifo.shift({});\n'.format(infmt, infmt)

            # write outputs
            oshape = thefdmt.hist_state_shape[iterno+1]

            write = '\n\n// Write outputs\n\n'
            for c in range(ncout):
                for d in range(ndout):
                    if (c, d) in sums_done: # if the sum was actually done - load it into the array
                        vout = fmt(iterno+1, c, d);
                    else:
                        vout = '0';  # set the few dangling outputs to zero

                    if ncout == 1:
                        write += '    out[{}] = {};\n'.format(d, vout)
                    else:
                        write += '    out[{}][{}] = {};\n'.format(d, c , vout)


            iters += iterstart + read  + do_sums + queuepush + write + '}\n\n'

        # sort queues by size
        sorted_queues = sorted(list(output_fifo_sizes.items()), key=lambda fsz: max(fsz[1]))
        nfifos = len(output_fifo_sizes)
        ngroups = int(np.ceil(float(nfifos)/float(fifos_per_group)))
        npad = ngroups*fifos_per_group - nfifos # number of FIFOS not to include in the first group
        assert 0 <= npad < fifos_per_group
        # NB: padding at the beginning rather than the end saves a bunch of memory
        group_sizes = np.zeros(ngroups, dtype=int)
        group_fifos = {}
        total_fifo_entries = 0
        num_lut16 = 0
        num_bram = 0
        word_size_bits = 16

        for fifo_enum, (fifo_name, fifo_sizes) in enumerate(sorted_queues):
            fifo_id = fifo_enum + npad # pad id from the beginning
            group_id = fifo_id // fifos_per_group
            group_offset = fifo_id % fifos_per_group
            fifo_size = max(fifo_sizes)
            total_fifo_entries += fifo_size
            if fifo_size >= compound_threshold:
                nlut = word_size_bits
            else:
                nlut = word_size_bits*((fifo_size + 16) // 16)

            if fifo_size >= compound_threshold:
                num_bram += (fifo_size + bram_depth -1 )//bram_depth
                
            print((fifo_name, fifo_size, nlut, num_bram))

            num_lut16 += nlut
            queuedecl += 'static FdmtFifo<{}, {}, {}> {}_fifo;\n'.format(fifo_size, group_id, group_offset, fifo_name);
            group_sizes[group_id] = max(group_sizes[group_id], fifo_size)
            gf = group_fifos.get(group_id, [])
            gf.append(fifo_name)
            group_fifos[group_id] = gf


        loadfun = ''
        header_loadfun_decl = ''

        cacheid = 0
        this_cache_depth = 0
        cache_sizes = []
        group_offsets = []
        group_cache_ids = []
        cache_groups = [[]]

        # load all groups functions
        load_all_fifos_decl = '''
void fdmt_load_fifos_from_cache(const group_cache_t input_cache[NUM_CACHES], group_cache_t output_cache[NUM_CACHES])'''
        header_loadfun_decl += load_all_fifos_decl + ';\n'

        loadallfun = load_all_fifos_decl + '''
{
#pragma HLS INLINE OFF
'''

        for group_id, fifos in group_fifos.items():
            # Work out cache sizes
            groupsz = group_sizes[group_id]
            if this_cache_depth + groupsz > max_cache_depth:
                cache_sizes.append(this_cache_depth)
                cacheid += 1
                this_cache_depth = 0
                cache_groups.append([])

            group_offsets.append(this_cache_depth)
            group_cache_ids.append(cacheid)
            cache_groups[-1].append(group_id)

            this_cache_depth += groupsz

            loadallfun += 'fdmt_load_fifos_group{group_id}(input_cache[FIFO_GROUP_CACHE_IDS[{group_id}]],output_cache[FIFO_GROUP_CACHE_IDS[{group_id}]]);\n'.format(**locals())

            # make load function
            loaddecl = '''void fdmt_load_fifos_group{group_id}(const group_cache_t input_cache,group_cache_t output_cache)'''.format(**locals())
            header_loadfun_decl += loaddecl + ';\n'
            loadfun += loaddecl + '''
    {{
        for(int i =0; i < FIFO_GROUP_NCLKS[{group_id}]; i++) {{
    #pragma HLS PIPELINE II=1

    '''.format(**locals())

            for fifo_name in fifos:
                loadfun += '        {fifo_name}_fifo.group_shift(i, input_cache, output_cache);\n'.format(**locals())

            loadfun += '''
        } // end for
    } // end fifo load

    '''
        loadallfun += '}\n\n'

        # Add final cache size
        cache_sizes.append(this_cache_depth)
        code = Code(preamble='', hfile=None)
        # Make cache load functions
        # make big cache load funciton
        load_all_func = code.cfunction('void fdmt_load_fifos_from_cache(const group_cache_t input_cache[NUM_CACHES], group_cache_t output_cache[NUM_CACHES])')
        load_all_func += '// Loads all caches in parallel'

        for cacheid, cache_groups in enumerate(cache_groups):
            cache_load_func = code.cfunction('void fdmt_load_cache{cacheid}(const group_cache_t input_cache, group_cache_t output_cache)'.format(**locals()))
            cache_load_func += '// Loads each group sequentially'
            load_all_func += 'fdmt_load_cache{cacheid}(input_cache[{cacheid}], output_cache[{cacheid}]);'.format(**locals())
            for groupid in cache_groups:
                cache_load_func += 'fdmt_load_fifos_group{groupid}(input_cache, output_cache);'.format(**locals())





        # output constants
        group_sizes_csep = ','.join(map(str, group_sizes))
        group_offsets_csep =  ','.join(map(str, group_offsets))
        cache_sizes_csep = ','.join(map(str, cache_sizes))
        group_cache_ids_csep = ','.join(map(str, group_cache_ids))
        num_caches = len(cache_sizes)
        total_nclks = sum(group_sizes)

        constants = '''
const int NC={f.n_f}; // Number of channels
const float FMIN = {f.f_min}; // Frequency of bottom channel (GHz)
const float FMAX = {f.f_max}; // Frequency of bottom channel (GHz)
const int ND = {f.max_dt}; // Number of output DM trials
const int ND_IN = {f.hist_state_shape[0][1]}; // number of input dm trials
const float DF = {f.d_f}; // channel interval (GHz)
const float BW = {f.bw}; // Total bandwidth (GHz)
'''.format(f=thefdmt)
        constants += '''
const int NGROUPS = {ngroups}; // total number of FIFO groups
const int NFIFOS_PER_GROUP = {fifos_per_group}; // Number of FIFOs in a group. Should = DRAM bus width = 64 in full system
const int MAX_CACHE_DEPTH = {max_cache_depth}; // Maximum depth of a cache block. Is the depth of the BRAM = 1024 in full system
const int FIFO_TOTAL_NCLKS = {total_nclks}; // Total number of CLKS require to clock in all the FIFOs serially
const int FIFO_GROUP_NCLKS[] = {{ {group_sizes_csep} }} ; // Number of clocks to load in each FIFO group
const int FIFO_GROUP_OFFSETS[] = {{ {group_offsets_csep} }}; // Offset address for each FIFO group within its cache block
const int NUM_CACHES = {num_caches}; // Total number of cache blocks
const int CACHE_SIZES[] = {{ {cache_sizes_csep} }}; // Depth of each individual cache (we can't always use all of a cache)
const int FIFO_GROUP_CACHE_IDS[] = {{ {group_cache_ids_csep} }}; // The ID of the cache each group will go in
const int TOTAL_FIFO_ENTRIES = {total_fifo_entries};
const int COMPOUND_THRESHOLD = {compound_threshold};
const int NUM_LUT16 = {num_lut16};
const int NUM_BRAM = {num_bram};
const int BRAM_DEPTH = {bram_depth};


typedef fdmt_t group_cache_t[MAX_CACHE_DEPTH][NFIFOS_PER_GROUP];
const char* const FDMT_NAME = "{self.root_file_name}";
'''.format(**locals())

        assert max(cache_sizes) <= max_cache_depth

        process_func_decl = 'void fdmt_process(fdmt_t in[ND_IN][NC], fdmt_t out[ND])'

        funcstart = process_func_decl + ''' {
        #pragma HLS PIPELINE II=16
        #pragma HLS INLINE off
        '''
        funcdecl = ''
        funcrun = ''
        lastiter = len(thefdmt.hist_nf_data) -1
        for iterno, nfd in enumerate(thefdmt.hist_nf_data):
            nextiter = iterno+1
            out_shape = thefdmt.hist_state_shape[nextiter]
            in_shape = thefdmt.hist_state_shape[iterno]
            if iterno != 0:
                funcdecl += '    fdmt_t d_iter{i}[{nd}][{nc}];\n'.format(i=iterno,nc=in_shape[0], nd=in_shape[1])
            if iterno == 0:
                funcrun += '    iteration{iterno}(in, d_iter{nextiter});\n'.format(**locals())
            elif iterno == lastiter:
                funcrun += '    iteration{iterno}(d_iter{iterno}, out);\n'.format(**locals())
            else:
                funcrun += '    iteration{iterno}(d_iter{iterno}, d_iter{nextiter});\n'.format(**locals())


        funcend = '\n\n}'
        fileend = "#endif"

        constants_com = ''
        for c in constants.split('\n'):
            constants_com += '// ' + c.strip() + '\n';

        if comment_constants:
            constants = constants_com

        cfile = '''
{preamble}
#include "fdmt.h"
{queuedecl}
{iters}
{funcstart}
{funcdecl}
{funcrun}
{funcend}

// FIFO loading funtions
{loadfun}
{code.code}
 '''.format(**locals())


        hfile = '''#ifndef _FDMT_PROCESS_H
#define _FDMT_PROCESS_H
{preamble}

// Constants
{constants}

// Iteration functions
{header_iterdecl}

// Processing functions
{process_func_decl};

// FIFO loading funtions
{code.header}


#endif
'''.format(**locals())

        self.cfile = cfile
        self.hfile = hfile
        self.preamble = preamble
        self.constants = constants
        self.queuedecl = queuedecl
        self.iters = iters
        self.funcstart = funcstart
        self.funcdecl = funcdecl
        self.funcrun = funcrun
        self.funcend = funcend
        self.loadfun = loadfun
        self.fileend = fileend

    @property
    def functional_file_name(self):
        '''Returns the file name that is only the functional stuff - not the cache depth or ff'''
        fn = 'fdmt_d{f.max_dt}_c{f.n_f}_f{f.f_min}'.format(f=self.thefdmt)
        return fn

    @property
    def root_file_name(self):
        root = 'fdmt_d{f.max_dt}_c{f.n_f}_f{f.f_min}_ff{self.fifos_per_group}_mcd{self.max_cache_depth}_iter3'.format(f=self.thefdmt, self=self)
        return root

    @property
    def header_file_name(self):
        return self.root_file_name + '.h'

    @property
    def code_file_name(self):
        return self.root_file_name + '.cpp'


    def write_files(self, target_dir='.'):
        print(('Writing to ', target_dir))
        hout = os.path.join(target_dir, self.header_file_name)
        with open(hout, 'w') as fout:
            print(("Writing to", hout))
            fout.write(self.hfile)

        cout = os.path.join(target_dir, self.code_file_name)
        with open(cout, 'w') as fout:
            print(("Writing to", cout));
            fout.write(self.cfile)

    def write_random_test_vectors(self, target_dir='.', seed=42):
        np.random.seed(seed)
        din = np.random.randint(-2**3, 2**3, size=(self.thefdmt.n_f, self.thefdmt.n_t)).astype(np.float32)
        dout = self.thefdmt(din)
        makedirs(target_dir, exist_ok=True)
        din.T.tofile(os.path.join(target_dir, self.functional_file_name+'.test.rand.in'))
        dout.T.tofile(os.path.join(target_dir, self.functional_file_name+'.test.rand.out'))

    def write_ones_test_vectors(self, target_dir='.'):
        din = np.ones((self.thefdmt.n_f, self.thefdmt.n_t), dtype=np.float32)
        dout = self.thefdmt(din)
        makedirs(target_dir, exist_ok=True)
        din.T.tofile(os.path.join(target_dir, self.functional_file_name+'.test.ones.in'))
        dout.T.tofile(os.path.join(target_dir, self.functional_file_name+'.test.ones.out'))

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Generates HLS code for an FDMT', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', action='store_true', help='Be verbose')
    parser.add_argument('--fch1', type=float, help='First channel (GHz)', default=1.12)
    parser.add_argument('--nchan', type=int, help='Number of channels', default=288)
    parser.add_argument('--chanbw', type=float, help='Channel bandwidtdh (MHz)', default=1.0)
    parser.add_argument('--nt', type=int, help='Time samples per block', default=256)
    parser.add_argument('--destdir', help='Desination directory for generated files', default='.')
    parser.add_argument('--test_dest_dir', help='Destination directory for test data', default='./testdata')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    # Typical center frequencies are 860 and 920 Mhz
    fc = 0.860 # center frequency GHz
    bw = 0.288 # bandwidth GHz
    Nd = 512 # number of DM trials
    Nchan= 256
    Nt = 256 # time block size
    Tint = 0.864e-3 # integration time - seconds
    Npol = 2 # input number of polarisations
    Npix = 256
    f1 = fc - bw/2.
    f2 = fc + bw/2.
    chanbw = 1e-3
    lam1 = constants.c/f1/1e9
    lam2 = constants.c/f2/1e9
    freqs = f1 + np.arange(Nchan)*chanbw
    lambdas = constants.c / (freqs*1e9)
    thefdmt = fdmt.Fdmt(f1, chanbw, Nchan, Nd, Nt)

    for nd, nchan in configs:
        thefdmt = fdmt.Fdmt(f1, chanbw, nchan, nd, Nt)
        dag = FdmtDag(thefdmt)

        if nd < 32:
            code = FdmtDagFileIter3(dag, fifos_per_group=nd, max_cache_depth=nd)
            code.write_files(values.destdir)
            code.write_ones_test_vectors(values.test_dest_dir)
            code.write_random_test_vectors(values.test_dest_dir)

        code = FdmtDagFileIter3(dag, fifos_per_group=32, max_cache_depth=1024)
        code.write_files(values.destdir)
        code.write_ones_test_vectors(values.test_dest_dir)
        code.write_random_test_vectors(values.test_dest_dir)


if __name__ == '__main__':
    _main()
