#!/usr/bin/env ParselTongue

################################################################################
# AIPS imports
################################################################################
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask
#from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData
#from AIPSData import AIPSUVData, AIPSImage, AIPSCat
#from AIPSTV import AIPSTV

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
parser.add_argument('-o', '--outfile', default="data", help="Output AIPS file")
parser.add_argument('FITSfile1', help="Input FITS data")
parser.add_argument('FITSfile2', help="Input FITS data")
parser.add_argument('FITSfile3', help="Input FITS data")
parser.add_argument('FITSfile4', help="Input FITS data")
args = parser.parse_args()

if args.user is not None: AIPS.userno = args.user

aipsDisk = 1

################################################################################
# Some useful functions
################################################################################
def fitld_uvfits(uvfitsfile, aipsdata):
    fitld = AIPSTask('fitld', version = aipsver)
    fitld.digicor = 1
    fitld.ncount = 1
    fitld.outdata = aipsdata
    if not os.path.exists(uvfitsfile):
        print uvfitsfile + " does not exist! Aborting."
        sys.exit()
    tempinfits = os.getcwd() + "/templink.uvfits"
    os.system("rm -f " + tempinfits)
    os.system("ln -s " + uvfitsfile + " " + tempinfits)
    fitld.datain = tempinfits
    fitld()
    os.system("rm -f " + tempinfits)

def gluUV(card_uvdata, glu_uvdata):
    vbglu = AIPSTask('vbglu', version = aipsver)

    if len(card_uvdata)!=4:
        print "Must supply 4 UVdata files to merge"
        exit()

    vbglu.indata = card_uvdata[0]
    vbglu.in2data = card_uvdata[1]
    vbglu.in3data = card_uvdata[2]
    vbglu.in4data = card_uvdata[3]

    vbglu.outdata = glu_uvdata

    vbglu()


def mergeIF(inuv, outname):    
    # Merge all IFs in input file to a single IF and reduce spectral resolution
    SpecAv = 24 # Amount of channels to average - should make dynamic
    
    IFcol = -1
    FreqCol = -1
    for i, t in enumerate(inuv.header['ctype']):
        if t=='IF':
            IFcol = i
        if t=='FREQ':
            FreqCol = i

    if IFcol==-1:
        print "Error: Could not find IF column in UV data"
        exit()
    if FreqCol==-1:
        print "Error: Could not find FREQ column in UV data"
        exit()

    num_if = inuv.header['naxis'][IFcol]
    num_freq = inuv.header['naxis'][FreqCol]

    if num_freq%SpecAv:
        print "Error: Cannot divide {} into {}. Aborting".format(SpecAv,num_freq)
        exit()

    avUV = AIPSUVData("TMPAV{}".format(pid), "UVDATA", 1,1)
    if avUV.exists(): avUV.zap()
        
    avspc = AIPSTask('avspc', version = aipsver)
    avspc.indata = inuv
    avspc.outdata = avUV
    avspc.avoption = 'SUBS'
    avspc.doacor = 1
    avspc.channel = SpecAv
    avspc()

    seq = 1
    while True:
        outData = AIPSUVData(outname, "UVDATA", aipsDisk, seq)
        if not outData.exists(): break
        seq += 1

    morif = AIPSTask('morif', version = aipsver)
    morif.indata = avUV
    morif.outdata = outData
    morif.npiece = 1/float(num_if)
    morif()

    avUV.zap()

    indxr = AIPSTask('indxr', version  = aipsver)
    indxr.indata = outData
    indxr()

    

    
################################################################################
# Main code
################################################################################

card_uvdata = []

fitsfile = [args.FITSfile1, args.FITSfile2, args.FITSfile3, args.FITSfile4]

for c, fits in enumerate(fitsfile):

    uvdata = AIPSUVData("TMPUV{}_{}".format(c, pid), "UVDATA", 1,1)
    if uvdata.exists(): uvdata.zap()

    #    inputfitsfile = "CRAFT_CARD{}.FITS".format(str(c))
    if not os.path.exists(fits):
        print "{} does not exist. Aborting".format(fits)
        sys.exit()
    
    # Load up the FITS file into AIPS
    fitld_uvfits(fits, uvdata)
    card_uvdata.append(uvdata)

gluUVdata = AIPSUVData("TMPGLU{}".format(pid), "UVDATA", 1,1)
if gluUVdata.exists(): glUVdata.zap()

gluUV(card_uvdata, gluUVdata)

mergeIF(gluUVdata, args.outfile)

for c in card_uvdata:
    c.zap()
gluUVdata.zap()

