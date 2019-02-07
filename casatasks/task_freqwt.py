from taskinit import *
import numpy as np
from scipy.interpolate import interp1d

def freqwt(vis, specfile, weightdata=True, cutoff=0.0):

    casalog.origin('freqwt')
    weights = np.loadtxt(specfile)
    if weights.ndim != 2 or weights.shape[1] != 2:
        casalog.post('Invalid weights shape {}'.format(weights.shape), priority='ERROR')
        return

    
    weightsf = weights[:, 0]
    weightsa = weights[:, 1]

    casalog.post('Weights from {}-{} GHz amplitude max/min/mean/std = {}/{}/{}/{}'.format(weightsf[0], weightsf[1], weightsa.max(), weightsa.min(), weightsa.mean(), weightsa.std()))
    interp = interp1d(weightsf, weightsa, kind='linear', fill_value='extrapolate')
    ms.open(vis)
    specwin = ms.getspectralwindowinfo()
    ms.close()
    if len(specwin) != 1:
        casalog.post('Cant handle more than 1 spectral window', priority='ERROR')
        return

    for w in specwin.values():
        chan_ghz = np.arange(w['NumChan'])*w['ChanWidth'] + w['Chan1Freq']
        chan_ghz /= 1e9
        w['ChanGhz'] = chan_ghz

        # interpolate onto data grid
        wamp = interp(chan_ghz)

        # set values less than threshold to 0
        wmask = wamp > cutoff
        wamp[~wmask] = 0

        # square and normalise weights spectrum
        wamp2 = wamp**2
        weightspec =  wamp2/sum(wamp2)

        # save in the dictionary
        w['weightspec'] = weightspec
        w['wamp'] = wamp
        w['wmask'] = wmask
        casalog.post('Got {} of {} valid weights above the cutoff'.format(sum(w['wmask']), len(wamp)))


    tb.open(vis, nomodify=False)
    for r in xrange(tb.nrows()):
        inw = tb.getcell('WEIGHT_SPECTRUM', r)
        # assume spectral window 0
        winno = '0'
        weightspec = specwin[winno]['weightspec']
        wamp = specwin[winno]['wamp']
        wmask = specwin[winno]['wmask']
        inw *= weightspec
        tb.putcell('WEIGHT_SPECTRUM', r, inw)

        if weightdata:
            ind = tb.getcell('DATA', r)
            ind[:, wmask] /= wamp[wmask]
            ind[:, ~wmask] = 0
            tb.putcell('DATA', r, ind)

    casalog.post('Weighted {} rows'.format(tb.nrows()))
                
    tb.addreadmeline('Modified weights from {}'.format(specfile))
    if weightdata:
        tb.addreadmeline('Modified data from {}'.format(specfile))

    tb.flush()
    tb.close()
    
    
    return True
