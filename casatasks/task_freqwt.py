from taskinit import *
import numpy as np
from scipy.interpolate import interp1d
import re


def load_spectrum(fname, polre=r'.*[_.]([XYIQUV]*)\.image-raster'):
    '''
    Loads spectrum from a text file - produced by CASA (somehow?)
    
    Assumes a file with lines starting with # as comments
    Each non-comment line contains freq in GHz and weight value separated by white space.
    Polarisation-dependant weights can be included in the file.
    Each polarisation should be prepended by a comment that satisfies the polarisation regexp. The first matched group of the regexp is used as the key for the return dictionaryy.
    If no polarisation header is in the file, 'I' is assumed
    :returns: A dictionary, keyed by polarisation type, a 2-D numpy array of (freqGHz, value)
    '''
    polre = re.compile(polre)
    currpol = 'I'
    out = {}
    out[currpol] = []
    last_freq = None
    with open(fname, 'rU') as f:
        for iline, line in enumerate(f):
            line = line.strip()
            if line.startswith('#'): # its a comment
                m = polre.match(line)
                if m is None: # no match - boring commment - continue to next line
                    continue
                else: # it's a new polarisation
                    currpol = m.groups(1)[0]
                    out[currpol] = []
                    last_freq = None
            else:
                bits = line.split()
                if len(bits) == 2:
                    freq, amp = map(float, bits)
                    if last_freq is not None and freq < last_freq:
                        s = 'Spectrum frequency is not monotonic on line {}! Please put only one spectrum in or split spectra by comment regex: {}'.format(iline+1, polre.pattern)
                        casalog.post(s, priority='ERROR')
                        raise ValueError(s)
                    last_freq = freq
                    out[currpol].append((freq, amp))

    # convert to numpy array
    out2 = {}
    for k, v in out.iteritems():
        out2[k] = np.array(v)
        assert out2[k].ndim == 2

    return out2
    
def plot_spec(spec):
    '''
    Plot the spectra loaded from load_spec
    '''
    from pylab import subplots, draw
    fig, ax = subplots(1,1)
    for pol, specd in spec.iteritems():
        f = specd[:, 0]
        a = specd[:, 1]
        ax.plot(f,a, label=pol)

    ax.set_xlabel('Frequency (GHz)')
    ax.legend()
    draw()

def freqwt(vis, specfile, weightdata=True, cutoff=0.0):

    import pylab
    casalog.origin('freqwt')
    try:
        spec = load_spectrum(specfile)
    except:
        casalog.post('Exception loading spectrum. No visibilities weighted', priority='ERROR')
        raise

    weight_interp = {}
    for pol, weightd in spec.iteritems(): #weights.iteritems():
        weightsf = weightd[:, 0]
        weightsa = weightd[:, 1]
        casalog.post('Weights file {} pol {} has {} frequencies from {:02f}-{:02f} GHz amplitude max/min/mean/std = {:02f}/{:02f}/{:02f}/{:02f}'
                     .format(specfile, pol, weightd.shape[0], weightsf[0], weightsf[-1], weightsa.max(), weightsa.min(), \
                             weightsa.mean(), weightsa.std()))
        weight_interp[pol] = interp1d(weightsf, weightsa, kind='linear', fill_value='extrapolate')

    ms.open(vis)
    specwin = ms.getspectralwindowinfo()
    ms.close()
    if len(specwin) != 1:
        casalog.post('Cant handle more than 1 spectral window', priority='ERROR')
        return

    for w in specwin.values():
        chan_ghz = (np.arange(w['NumChan'])*w['ChanWidth'] + w['Chan1Freq'])/1e9
        w['ChanGhz'] = chan_ghz

        # interpolate onto data grid
        interp = weight_interp['I'] # only support stokes I for now
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
        casalog.post('Got {} of {} valid weights above the cutoff of {}'.format(sum(w['wmask']), len(wamp), cutoff))


    tb.open(vis, nomodify=False)
    nblank = 0
    for r in xrange(tb.nrows()):
        inw = tb.getcell('WEIGHT_SPECTRUM', r)
        # assume spectral window 0
        winno = '0'
        weightspec = specwin[winno]['weightspec']
        if np.all(weightspec.flatten() == 0):
            weightspec += 1
            nblank += 1
        wamp = specwin[winno]['wamp']
        wmask = specwin[winno]['wmask']
        inw *= weightspec
        tb.putcell('WEIGHT_SPECTRUM', r, inw)

        if weightdata:
            ind = tb.getcell('DATA', r)
            # ind shape is (pol, freq)
            ind[:, wmask] /= wamp[wmask]
            ind[:, ~wmask] = 0
            tb.putcell('DATA', r, ind)

    casalog.post('Weighted {} rows'.format(tb.nrows()))
                
    tb.addreadmeline('Modified weights from {}. {} blank rows '.format(specfile, nblank))
    casalog.post('Replaced weights in {} rows'.format(nblank))
    if weightdata:
        tb.addreadmeline('Modified data from {}'.format(specfile))

    tb.flush()
    tb.close()
    
    
    return True
