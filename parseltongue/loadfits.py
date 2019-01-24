#!/usr/bin/env ParselTongue

################################################################################
# AIPS imports
################################################################################
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
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
    aipsver = os.environ['AIPSVER']
except KeyError:
    aipsver = '31DEC18'
AIPS.userno = 702

parser = argparse.ArgumentParser()
parser.add_argument('-u', '--user', help="AIPS user number", type=int)
parser.add_argument('-D', '--disk', default=1, help="AIPS disk", type=int)
parser.add_argument('-o', '--outfile', default="data", help="Output AIPS file")
parser.add_argument('-f', '--fitsfileoutname', help="Output FITS file, default doesn't write")
parser.add_argument('-a', '--antlist', help="Force antenna list")
parser.add_argument('-s', '--specav', default=18, help="Spectral Averaging", type=int)
parser.add_argument("-m", "--nomerge", default=False, action="store_true", help="Don't Merge IFs")
parser.add_argument("-d", "--delete", default=False, action="store_true", help="Delete AIPS file if writing FITS out")
parser.add_argument('fitsfile', nargs='+', help="Input FITS data")
args = parser.parse_args()

if args.user is not None: AIPS.userno = args.user

aipsDisk = args.disk

fitsfileoutname = args.fitsfileoutname
if fitsfileoutname is not None:
    if not fitsfileoutname[0] == '/':
        fitsfileoutname = os.getcwd() + '/' + args.fitsfileoutname
        if os.path.exists(fitsfileoutname):
        print "Warning: {} already exists. Aborting".format(fitsfileoutname)
        sys.exit()

################################################################################
# Some useful functions
################################################################################
def fitld_uvfits(uvfitsfile, aipsdata, antlist=None, specAv=None):
    # Load FITS file and immediatetly reduce spectral resolution
    if specAv is None: specAv = 18 # Amount of channels to average

    if not os.path.exists(uvfitsfile):
        print uvfitsfile + " does not exist! Aborting."
        sys.exit()

    if specAv==0:
        loadUV = aipsdata
    else:
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

    if antlist is not None:
        antname = [x.strip() for x in antlist.split(',')]
        fitld.antname = AIPSList(antname)
    
    fitld()
    os.system("rm -f " + tempinfits)

    if specAv>0:
        FreqCol = None
        for i, t in enumerate(loadUV.header['ctype']):
            if t=='FREQ':
                FreqCol = i
                break

        if FreqCol is None:
            print "Error: Could not find FREQ column in UV data"
            exit()

        num_freq = loadUV.header['naxis'][FreqCol]
        if num_freq%specAv:
            print "Error: Cannot divide {} into {}. Aborting".format(specAv,num_freq)
            exit()

        avspc = AIPSTask('avspc', version = aipsver)
        avspc.indata = loadUV
        avspc.outdata = aipsdata
        avspc.avoption = 'SUBS'
        avspc.doacor = 1
        avspc.channel = specAv
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
    ntmp = 1
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
            thisout = AIPSUVData("GLU{}-{}".format(pid, ntmp), "UVDATA", aipsDisk, 1)
            if thisout.exists(): thisout.zap()
            vbglu.outdata = thisout
            ntmp += 1

        vbglu()
        if lastMerge is not None:
            lastMerge.zap()
        lastMerge = thisout

    if thisout is not None:
        thisout.zap()

        
def mergeIF(inuv, outData):    
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

single = len(args.fitsfile)==1

seq = 1
while True:
    outdata = AIPSUVData(args.outfile, "UVDATA", aipsDisk, seq)
    if not outdata.exists(): break
    seq += 1


for c, fits in enumerate(args.fitsfile):

    if single and args.nomerge:
        uvdata = outdata
    else:
        uvdata = AIPSUVData("TMPUV{}_{}".format(c+1, pid), "UVDATA", aipsDisk, 1)
        if uvdata.exists(): uvdata.zap()

    if not os.path.exists(fits):
        print "{} does not exist. Aborting".format(fits)
        sys.exit()
    
    # Load up the FITS file into AIPS
    fitld_uvfits(fits, uvdata, args.antlist, args.specav)
    card_uvdata.append(uvdata)

if single:
    gluUVdata = card_uvdata[0]    
else:
    gluUVdata = AIPSUVData("TMPGLU{}".format(pid), "UVDATA", aipsDisk, 1)
    if gluUVdata.exists(): glUVdata.zap()
    gluUV(card_uvdata, gluUVdata)

if not args.nomerge:
    mergeIF(gluUVdata, outdata)
    for c in card_uvdata:
        c.zap()
    if not single: gluUVdata.zap()
else:
    if not single:
        for c in card_uvdata:
            c.zap()
    
if args.fitsfileoutname is not None:
    if not args.fitsfileoutname[0] == '/':
        args.fitsfileoutname = os.getcwd() + '/' + args.fitsfileoutname
    fittp = AIPSTask("fittp")
    fittp.indata = outdata
    fittp.dataout = args.fitsfileoutname
    fittp()
    if args.delete:
        outdata.zap()
