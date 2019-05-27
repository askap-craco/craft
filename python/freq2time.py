#!/usr/bin/env python
"""
reconstructs FFFF via de-rippling, coherent de-dispersing, and IFFT-ing.
"""

import numpy as np
from scipy.interpolate import interp1d
import os
from scipy import io

def generate_deripple(nfft,res):
    N = 1536
    OS_De = 27.
    OS_Nu = 32.
    dr = './ADE_R6_OSFIR.mat'
    h = io.loadmat(dr)['c'][0]
    passbandLength = int(((nfft / 2) * OS_De) / OS_Nu)
    multiple = int(float(N)/res)
    deripple=abs(np.fft.fft(h,multiple*passbandLength*2))
    return deripple

def deripple(FFFF, fftLength = 1048576, quiet=False):
    if not quiet:
        print('derippling....')
    FFFF= FFFF[0,:,0]

    # ASKAP Parameters
    N = 1536
    OS_De = 27.
    OS_Nu = 32.
    passbandLength = int(((fftLength / 2) * OS_De) / OS_Nu)

    # de-ripple coefficients
    res = 6 # resolution of the deripple coefficients
    temp = generate_deripple(fftLength,res)
    interp = interp1d(res*np.arange(len(temp)),temp)
    deripple = np.ones(passbandLength+1)/abs(interp(np.arange(passbandLength+1)))

    for chan in range(336):
    #print(chan)
        for ii in range(passbandLength):
            FFFF[ii+chan*passbandLength*2] = FFFF[ii+chan*passbandLength*2]*deripple[passbandLength-ii]
            FFFF[passbandLength+ii+chan*passbandLength*2] = FFFF[passbandLength+ii+chan*passbandLength*2]*deripple[ii]
    
    return FFFF

def coh_dedisp(FFFF, DM, f_mid=1320.5, bw=336, quiet=False):
    nSam = len(FFFF)

    # ASKAP Parameters
    f_start = f_mid - float(bw)/2 #1153.
    f_stop = f_mid + float(bw)/2 #1488.
    
    if not quiet:
        print('dedispersing....')
    dedisperse = np.exp(2j*np.pi*DM/2.41e-4*np.array([(f-f_mid)**2/f_mid**2/f*1e6 for f in np.linspace(f_stop,f_start,nSam)]))
    #print('dedispersing wrong....')
    #dedisperse = np.exp(2j*np.pi*DM*4150*np.array([(1/f**2-1/f_mid**2)*f*1e6 for f in np.linspace(f_stop,f_start,nSam)]))

    FFFF *= dedisperse

    return FFFF

def ifft_long(FFFF,quiet=False):
    if not quiet:
        print('ifft-ing....')
    t_series= np.fft.ifft(np.fft.fftshift(FFFF))
    return t_series

def reconstruct(fn, fftLength, DM, f0 = 1320.5, bw=336, quiet=False):
    FFFF = np.load(fn)
    FFFF = deripple(FFFF,fftLength,quiet)
    FFFF = coh_dedisp(FFFF,DM,f0,bw,quiet)
    t_series = ifft_long(FFFF,quiet)
    return t_series


def _main():
    import time
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    t0 = time.time()
    parser = ArgumentParser(description='Script description', formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--fn', help="FFFF file directory", default=None)
    parser.add_argument('-d', '--DM', type=float, help='Dispersion measure', default=None)
    parser.add_argument('--f0', type=float, help='Central frequency', default=1320.5)
    parser.add_argument('-o', '--outfile', help="Output time series file directory")
    parser.add_argument('-l', '--fftlength', type=int, help="FFT length", default=1048576)
    parser.add_argument('--no_dr', help="Don't deripple", action='store_true', default=False)
    parser.add_argument('--no_dd', help="Don't dedisperse", action='store_true', default=False)
    parser.add_argument('--no_ifft', help="Don't ifft", action='store_true', default=False)
    parser.add_argument('-q', help="Quiet mode", action='store_true', default=False)
    values = parser.parse_args()

    if values.no_dr or values.no_dd or values.no_ifft:
        FFFF = np.load(values.fn)
        if values.no_dr and values.no_dd:
            print('No derippling, no dedispersing.')
            FFFF= FFFF[0,:,0]
        elif values.no_dr:
            print('No derippling')
            FFFF= FFFF[0,:,0]
            FFFF = coh_dedisp(FFFF,values.DM,f_mid=values.f0, quiet=values.q)
        else:
            print('No dedispersing')
            FFFF = deripple(FFFF,values.fftlength, quiet=values.q)
        if values.no_ifft:
            print('No ifft. Saving FFFF')
            t_series = FFFF
        else:
            t_series = ifft_long(FFFF,quiet=values.q)
    else:
        print('Reconstructing...')
        t_series = reconstruct(values.fn, values.fftlength, values.DM, f0=values.f0, quiet=values.q)
    
    if values.outfile is not None:
        print('output saved to: '+values.outfile)
        np.save(values.outfile,t_series)
    print('freq2time.py running time: '+str(time.time()-t0))
        
if __name__ == '__main__':
    _main()
    
