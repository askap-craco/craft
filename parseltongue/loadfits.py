#!/usr/bin/env ParselTongue

################################################################################
# AIPS imports
################################################################################
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask
from AIPSData import AIPSUVData

################################################################################
# General python imports
################################################################################
import argparse, sys, os, math

pid = os.getpid()

################################################################################
# Set up AIPS stuff
################################################################################
try:
    aipsver = os.environ['PSRPIAIPSVER']
except KeyError:
    aipsver = '31DEC18'
AIPS.userno = 702

parser = argparse.ArgumentParser()
parser.add_argument('-u', '--user', help="AIPS user number", type=int)
parser.add_argument('-d', '--disk', default=1, help="AIPS disk", type=int)
parser.add_argument('-o', '--outfile', default="data", help="Output AIPS file")
parser.add_argument('fitsfile', nargs='+', help="Input FITS data")
args = parser.parse_args()

if args.user is not None: AIPS.userno = args.user

aipsDisk = args.disk

################################################################################
# Some useful functions
################################################################################
def fitld_uvfits(uvfitsfile, aipsdata):
    # Load FITS file and immediatetly reduce spectral resolution
    SpecAv = 18 # Amount of channels to average - should make dynamic

    if not os.path.exists(uvfitsfile):
        print uvfitsfile + " does not exist! Aborting."
        sys.exit()

    loadUV = AIPSUVData("TMPLOAD{}".format(pid), "UVDATA", aipsDisk, 1)
    if loadUV.exists(): loadUV.zap()
        
    fitld = AIPSTask('fitld', version = aipsver)
    fitld.digicor = 1
    fitld.ncount = 1
    fitld.outdata = loadUV
    tempinfits = os.getcwd() + "/templink.uvfits"
    os.system("rm -f " + tempinfits)
    os.system("ln -s " + uvfitsfile + " " + tempinfits)
    fitld.datain = tempinfits
    fitld()
    os.system("rm -f " + tempinfits)

    FreqCol = None
    for i, t in enumerate(loadUV.header['ctype']):
        if t=='FREQ':
            FreqCol = i
            break

    if FreqCol is None:
        print "Error: Could not find FREQ column in UV data"
        exit()

    num_freq = loadUV.header['naxis'][FreqCol]
    if num_freq%SpecAv:
        print "Error: Cannot divide {} into {}. Aborting".format(SpecAv,num_freq)
        exit()

    avspc = AIPSTask('avspc', version = aipsver)
    avspc.indata = loadUV
    avspc.outdata = aipsdata
    avspc.avoption = 'SUBS'
    avspc.doacor = 1
    avspc.channel = SpecAv
    avspc()

    loadUV.zap()

def gluUV(card_uvdata, glu_uvdata):
    vbglu = AIPSTask('vbglu', version = aipsver)

    nfits =  len(card_uvdata)
    if nfits<2:
        print "Glu'ing a single datfile makes no sense. Aborting"
        sys.exit()

    lastMerge = None
    thisout = None
    nleft = nfits
    while nleft>0:
        if lastMerge != None:
            vbglu.indata = lastMerge
        else:
            vbglu.indata = card_uvdata[nfits-nleft]
            nleft -= 1
        vbglu.in2data = card_uvdata[nfits-nleft] # Will always be set
        nleft -= 1
        if nleft>0:
            vbglu.in3data = card_uvdata[nfits-nleft]
            nleft -= 1
            if nleft>0:
                vbglu.in4data = card_uvdata[nfits-nleft]
                nleft -= 1
            else:
                vbglu.in4name = ""
                vbglu.in4class = ""
                vbglu.in4seq = 0
                vbglu.in4disk = 0
        else:
            vbglu.in3name = ""
            vbglu.in3class = ""
            vbglu.in3seq = 0
            vbglu.in3disk = 0
            vbglu.in4name = ""
            vbglu.in4class = ""
            vbglu.in4seq = 0
            vbglu.in4disk = 0
        
        if nleft==0:
            vbglu.outdata = glu_uvdata
            thisout = None
        else:
            thisout = AIPSUVData("GLU{}".format(pid), "UVDATA", aipsDisk, 1)
            if thisout.exists(): thisout.zap()
            vbglu.outdata = thisout

        vbglu()
        if nleft>0:
            if lastMerge is not None:
                lastMerge.zap()
                pass
            lastMerge = thisout
    if thisout is not None:
        thisout.zap()
def mergeIF(inuv, outname):    
    # Merge all IFs in input file to a single IF and reduce spectral resolution
    IFcol = -1
    for i, t in enumerate(inuv.header['ctype']):
        if t=='IF':
            IFcol = i
            break
        
    if IFcol==-1:
        print "Error: Could not find IF column in UV data"
        exit()

    num_if = inuv.header['naxis'][IFcol]

    seq = 1
    while True:
        outData = AIPSUVData(outname, "UVDATA", aipsDisk, seq)
        if not outData.exists(): break
        seq += 1

    morif = AIPSTask('morif', version = aipsver)
    morif.indata = inuv
    morif.outdata = outData
    morif.npiece = 1/float(num_if)
    morif()

    indxr = AIPSTask('indxr', version  = aipsver)
    indxr.indata = outData
    indxr()
    
################################################################################
# Main code
################################################################################

card_uvdata = []

for c, fits in enumerate(args.fitsfile):

    uvdata = AIPSUVData("TMPUV{}_{}".format(c+1, pid), "UVDATA", aipsDisk, 1)
    if uvdata.exists(): uvdata.zap()

    if not os.path.exists(fits):
        print "{} does not exist. Aborting".format(fits)
        sys.exit()
    
    # Load up the FITS file into AIPS
    fitld_uvfits(fits, uvdata)
    card_uvdata.append(uvdata)

gluUVdata = AIPSUVData("TMPGLU{}".format(pid), "UVDATA", aipsDisk, 1)
if gluUVdata.exists(): glUVdata.zap()

gluUV(card_uvdata, gluUVdata)

mergeIF(gluUVdata, args.outfile)

for c in card_uvdata:
    c.zap()
gluUVdata.zap()

