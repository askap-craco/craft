#!/usr/bin/env python
"""
Template for making scripts to run from the command line

Copyright (C) CSIRO 2017
"""
import pylab
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
import sys
import logging
import h5py
import glob

__author__ = "Keith Bannister <keith.bannister@csiro.au>"

class OdcFile(object):
    def __init__(self, fname, intoff=1):
        self.fname = fname
        self.f = h5py.File(self.fname, 'r')
        cdata = self.f['CCdata'].value
        (self.ntimes, self.nsubbands, self.nports) = cdata.shape
        self.intoff = intoff
        self.ntimes -= intoff
        freqs = self.f['skyFrequency'] # (nsubbands, ntimes)
        assert freqs.shape == (self.ntimes + self.intoff,  self.nsubbands)
        # start at the second itegration, the first is often junk
        chan0vtime = list(freqs[self.intoff:,0])
        # number of times per integration = how many times until we see the first channel again
        self.ntimes_per_integration = chan0vtime[1:].index(chan0vtime[0]) + 1
        self.nint = self.ntimes/self.ntimes_per_integration
        self.nfreq = self.ntimes_per_integration*self.nsubbands
        assert 0 < self.ntimes_per_integration < self.ntimes

    def get_trange(self, intno):
        assert 0 <= intno < self.nint
        trange = slice(intno*self.ntimes_per_integration + self.intoff,(intno+1)*self.ntimes_per_integration + self.intoff)
        return trange


    def data(self, intno):
        ''' 
        Returns the data and frequecy values for integration numeber i
        '''
        trange = self.get_trange(intno)
        nfreq = self.nfreq
        nports = self.nports
        f = self.f

        count = f['CCcount'][trange, :].reshape(nfreq)
        data = f['CCdata'][trange,:,:].reshape(nfreq,nports)
        status = f['CCstatus'][trange, :].reshape(nfreq)
        # todo: check status
        freq = f['skyFrequency'][trange,:].reshape(nfreq)
        idx = np.argsort(freq)
        mask = status[idx] == 0

        freq = freq[idx]
        data = data[idx, :]/count[idx, None]
        data = np.ma.masked_array(data, mask=np.tile(mask, (nports, 1)))

        return freq, data
        
    def metadata(self, intno):
        '''
        Returns average metadata for an integration number in a dictionary
        '''

        f = self.f
        trange = self.get_trange(intno)
        d = {}
        for s in ('azimuth','elevation','decJ2000','raJ2000','rollAngle','bat'):
            d[s] = f[s].value[trange].mean()

        return d
        
    @property
    def antenna(self):
        return self.f.attrs['antenna']

    def close(self):
        try:
            self.f.close()
        except:
            pass

    def __del__(self):
        self.close()

def group_by_ant(files):
    ref_files = {}
    for ref_file in files:
        try:
            refin = OdcFile(ref_file)
            ant = refin.antenna
            all_refs = ref_files.get(ant, [])
            all_refs.append(refin)
            ref_files[ant] = all_refs
            logging.debug('Opened reference files %s for antenna %s which has %d reference files', ref_file, ant, len(all_refs))
        except:
            logging.exception('Couldnt open file %s', ref_file)
            

    return ref_files


npaf = 188

def _main():
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Be verbose')
    parser.add_argument('-r', '--refdir', help='Reference directory')
    parser.add_argument('-s', '--show', action='store_true', help='Show plots')


    parser.add_argument(dest='files', nargs='+')
    parser.set_defaults(verbose=False)
    values = parser.parse_args()
    if values.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    # parse the reference direcotry:
    ref_files = group_by_ant(glob.glob(os.path.join(values.refdir, '*.cal.hdf5')))
    logging.debug('Got reference files for antennas %s', str(list(ref_files.keys())))

    ant_files = group_by_ant(values.files)


    for ant, ant_files in ant_files.items():
        logging.debug('Processing files %s', ant_files)
        ref_odc = ref_files.get(ant, [])
        if len(ref_odc) > 0:
            ref_odc = ref_odc[0]
            summarise_ref(ant_files, ref_odc, values)
        else:
            logging.info('No refrerence for antenna: %s. Skipping' % ant)

def summarise_ref(ant_files, ref_odc, values):
    #ref_freq, ref_data = odcf.data(0)

    normd = []
    all_metadata = []
    chanrange = slice(None,None)
    fild = []
    delayfits = {}

    for odcf in ant_files:

        f = odcf.f
        nint = odcf.nint
        ref_freq, ref_data = ref_odc.data(0)
        ref_freq = ref_freq[chanrange]
        ref_data = ref_data[chanrange]
        print('Running', odcf.fname)


        for i in range(nint):
            freq, data = odcf.data(i)
            if np.any(freq != ref_freq):
                #print 'Bleach freqs to match!', (freq - ref_freq)
                print(freq)
                print(ref_freq)
                fig,ax = pylab.subplots(1,1)
                mask = freq != ref_freq
                ax.plot(freq, ref_freq)
                ax.plot(freq[mask], ref_freq[mask], 'ro')
                ax.set_xlabel('Current frequency')
                ax.set_ylabel('Reference freqency')
                pylab.show()


            freq = freq[chanrange]
            data = data[chanrange]
            metadata = odcf.metadata(i)
            all_metadata.append(metadata)

            #assert np.all(freq == ref_freq)

            norm_data = data/ref_data
            normd.append(norm_data)
            dofilter = False
            if dofilter:
                dnorm4 = norm_data.copy()
                for prt in range(npaf):
                    d = data[:, prt]
                    delayfit = delayfits.get(prt, None)
                    if delayfit is None:
                        delayfit = np.polyfit(freq, np.unwrap(np.angle(d)), 1)
                        delayfits[prt] = delayfit
                        
                        d3 = abs(d)*np.exp(1j*(np.angle(d) - np.polyval(delayfit, freq)))
                    d3f = np.fft.fft(d3)
                    d3f[3:] = 0
                    dnorm4[:, prt] = np.fft.ifft(d3f)
                
                    fild.append(dnorm4)

            if values.show:
                pts = np.arange(188)
                
                fig, ax = pylab.subplots(1,2)
                ax[0].plot(freq, abs(data))
                ax[1].plot(freq, abs(ref_data))
                
                fig, ax = pylab.subplots(1,3)
                fd = np.fft.fft(abs(data[:, :]), axis=0)
                fd[0, :] = 0
                
                ax[0].plot(freq, abs(norm_data[:,pts ]))
                ax[1].plot(freq, np.degrees(np.angle(norm_data[:, pts])))
                ax[0].set_ylabel('Abs(norm)')
                ax[1].set_ylabel('angle(norm)')
                ax[2].semilogy(np.fft.fftshift(abs(fd)))
                
                pylab.figure()
                pylab.imshow(np.log10(np.fft.fftshift(abs(fd), axes=0)), aspect='auto')
                
                pylab.figure()
                pylab.imshow(abs(data), vmin=0.5,vmax=1.5, aspect='auto')
                pylab.xlabel('Port number')
                pylab.ylabel('Channel')
                pylab.title('ODC amplitude normalised by port 48')
                pylab.show()

    fild = np.array(fild)
    normd = np.array(normd)
    #normd = fild
    normd /= normd[-1, :, :]

    print('Fild', fild.shape, 'normd', normd.shape)
    
    bats = np.array([m['bat'] for m in all_metadata])
    times = (bats - bats[0])/1e6/60
    refpt = 42
    refchan = 38
    fig, ax = pylab.subplots(1,2, sharex=True)
    ax[0].plot(freq, abs(normd[:,:,refpt]).T)
    ax[0].set_title('Amplitude change port {}'.format(refpt+1))
    ax[1].plot(freq, np.degrees(np.angle(normd[:,:,refpt])).T)
    ax[1].set_title('Phase change (deg) port {}'.format(refpt+1))
    ax[1].set_xlabel('Freq (MHz)')
    fout = odcf.fname + '.ampphvf_p{}.png'.format(refpt+1)
    fig.suptitle(fout)
    pylab.savefig(fout)



    fig, ax = pylab.subplots(1,2, sharex=True)
    ax[0].imshow(abs(normd[:,:,refpt]).T, aspect='auto')
    ax[0].set_title('Amplitude change port {}'.format(refpt+1))
    ax[1].imshow(np.degrees(np.angle(normd[:,:,refpt])).T,  aspect='auto')
    ax[1].set_title('Phase change (deg) port {}'.format(refpt+1))
    ax[0].set_ylabel('Channel')
    ax[0].set_xlabel('Integration')
    fout = odcf.fname + '.ampphim_p{}.png'.format(refpt+1)
    fig.suptitle(fout)
    pylab.savefig(fout)

    fig, ax = pylab.subplots(1,2, sharex=True)

    ax[0].plot(times, abs(normd[:,:,refpt]))
    ax[0].set_title('Amplitude change port {}'.format(refpt+1))
    ax[1].plot(times, np.degrees(np.angle(normd[:,:,refpt])))
    ax[1].set_title('Phase change (deg) port {}'.format(refpt+1))
    ax[1].set_xlabel('Time(min)')
    ax[0].set_xlabel('Time(min)')

    fout = odcf.fname + '.ampphvt_p{}.png'.format(refpt+1)
    fig.suptitle(fout)
    pylab.savefig(fout)

    fig, ax = pylab.subplots(1,2, sharex=True)
    ax[0].plot(times, abs(normd[:,refchan,0:npaf]))
    ax[0].set_title('Amplitude change channel {}'.format(refchan))
    ax[1].plot(times, np.degrees(np.angle(normd[:,refchan,0:npaf])))
    ax[1].set_title('Phase change (deg) channel {}'.format(refchan))
    ax[0].set_xlabel('Time(min)')
    ax[1].set_xlabel('Time(min)')

    fout = odcf.fname + '.amphvt.c{}.png'.format(refchan)
    fig.suptitle(fout)
    pylab.savefig(fout)

    fig, ax = pylab.subplots(1,2)
    ax[0].imshow(abs(normd[:,refchan,0:npaf].T), aspect='auto')
    ax[0].set_title('Amplitude change channel {}'.format(refchan))
    ax[1].imshow(np.degrees(np.angle(normd[:,refchan,0:npaf].T)), aspect='auto')
    ax[1].set_title('Phase change (deg) channel {}'.format(refchan))
    ax[0].set_xlabel('Time(min)')
    ax[1].set_xlabel('Time(min)')
    ax[0].set_ylabel('Port')

    fig, ax = pylab.subplots(1,2)

    ax[0].imshow(abs(normd[:,:,0:npaf].std(axis=0)).T, aspect='auto')
    ax[0].set_title('std of time')
    ax[1].imshow(np.degrees(np.angle(normd[:,:,0:npaf].std(axis=0))).T, aspect='auto')
    ax[1].set_title('std of time')
    ax[0].set_xlabel('Freq')
    ax[1].set_xlabel('Freq')
    ax[0].set_ylabel('Port')



    fig, ax = pylab.subplots(1,2)

    ax[0].imshow(abs(normd[:,:,0:npaf].mean(axis=2)), aspect='auto')
    ax[0].set_title('Mean of ports')
    ax[1].imshow(np.degrees(np.angle(normd[:,:,0:npaf])).mean(axis=2), aspect='auto')
    ax[1].set_title('Phase meanof ports')
    ax[0].set_xlabel('Freq')
    ax[1].set_xlabel('Freq')
    ax[0].set_ylabel('time')


    fig, ax = pylab.subplots(1,2)

    ax[0].plot(abs(normd[:,:,0:npaf].mean(axis=2)).T)
    ax[0].set_title('Mean of ports')
    ax[1].plot(np.degrees(np.angle(normd[:,:,0:npaf])).mean(axis=2).T)
    ax[1].set_title('Phase meanof ports')
    ax[0].set_xlabel('Freq')
    ax[1].set_xlabel('Freq')
    ax[0].set_ylabel('time')


    


    pylab.show()
                   

    

if __name__ == '__main__':
    _main()
