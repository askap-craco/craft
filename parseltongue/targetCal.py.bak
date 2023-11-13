#!/usr/bin/env ParselTongue
#Imports ########################################################
from AIPS import AIPS
from AIPSTask import AIPSTask
from AIPSData import AIPSUVData
import argparse, os, sys, glob

################################################################################
# Global variables and option parsing
################################################################################
try:
    aipsver = os.environ['PSRVLBAIPSVER']
except KeyError:
    aipsver = '31DEC18'

parser = argparse.ArgumentParser()
parser.add_argument("-r", "--refant", type=int, default=3,
                  help="The reference antenna to use.")
parser.add_argument("-u", "--userno", type=int, default=2,
                  help="The AIPS user number")
parser.add_argument("-s", "--sourcename", default="CRAFTSRC",
                  help="The name of the source in the FITS files")
parser.add_argument("-F", "--flagfile", 
                    help="Flag file to apply to target data, if desired.")
parser.add_argument('aipscalfile', help="AIPS calibrator file")
parser.add_argument('aipstargetfile', help="AIPS  target file")

args = parser.parse_args()

AIPS.userno     = args.userno
refant          = args.refant
aipsdisk        = 1

(name, aipsclass, seq) = args.aipscalfile.split('.')
caldata = AIPSUVData(name, aipsclass, aipsdisk, int(seq))
if not caldata.exists():
    print args.aipscalfile, "does not exists! Aborting"
    sys.exit()

(name, aipsclass, seq) = args.aipstargetfile.split('.')
targetdata = AIPSUVData(name, aipsclass, aipsdisk, int(seq))
if not targetdata.exists():
    print args.aipstargetfile, "does not exists! Aborting"
    sys.exit()
    
################# AIPS TASKS ################

def flag_data(uvdataset, filename, outflagver=1):
    uvflg = AIPSTask('uvflg', version = aipsver)
    uvflg.indata = uvdataset
    uvflg.intext = filename
    uvflg.outfgver = outflagver
    uvflg()

def applySN(uvdataset, snversion, clversion, refant=1, interpol='SELN'):
    clcal = AIPSTask('clcal', version = aipsver)
    clcal.indata = uvdataset
    clcal.sources[1:] = clcal.soucode = clcal.calsour[1:] = ''
    clcal.opcode = 'CALI'
    clcal.interpol = interpol
    clcal.doblank = 1
    clcal.dobtween = -1
    clcal.snver = snversion
    clcal.gainver = clversion
    clcal.refant=refant
    clcal()

def copyTab(indataset, tabletype, inver, outdataset, outver):
    tacop = AIPSTask('tacop', version = aipsver)
    tacop.indata = indataset
    tacop.inext = tabletype
    tacop.ncount = 1
    tacop.invers = inver
    tacop.outdata = outdataset
    tacop.outvers = outver
    tacop()
    

# Flag the calibrator data, if desired
if args.flagfile is not None:
    flag_data(targetdata, args.flagfile, 1)


# Copy two SN tables and BP table (shoud check these exists)

copyTab(caldata, "SN", 1, targetdata, 1)
copyTab(caldata, "SN", 2, targetdata, 2)
copyTab(caldata, "BP", 1, targetdata, 1)

# Apply SN tables 
applySN(targetdata, 1, 1, refant)
applySN(targetdata, 2, 2, refant)
