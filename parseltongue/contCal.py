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
parser.add_argument("-r", "--refant", type=int, default=1,
                  help="The reference antenna to use.")
parser.add_argument("-u", "--userno", type=int, default=2,
                  help="The AIPS user number")
parser.add_argument("-s", "--sourcename", default="CRAFTSRC",
                  help="The name of the source in the FITS files")
parser.add_argument("-f", "--flux", type=float, default=9.5, # 0407 flux
                  help="Calibrator flux in Jy,  Defaulted to correct value for 0407")
parser.add_argument("-F", "--flagfile", 
                  help="Flag file to apply to calibrator data only, if desired. Used to ensure RFI doesn't corrupt FRING or BPASS.")
parser.add_argument('aipsfile', help="AIPS  file ")

args = parser.parse_args()

AIPS.userno     = args.userno
refant          = args.refant
aipsdisk        = 1

(name, aipsclass, seq) = args.aipsfile.split('.')

caldata = AIPSUVData(name, aipsclass, aipsdisk, int(seq))

if not caldata.exists():
    print args.aipsfile, "does not exists! Aborting"
    sys.exit()

################# AIPS TASKS ################

def run_fring(aipsdata, refant=1, solint=1):
    fring = AIPSTask('fring', version = aipsver)
    fring.indata = aipsdata
    fring.refant = refant
    fring.solint = solint
    fring.dparm[9] = 1
    fring()

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


def run_bpass(uvdataset, srcname, clversion, refant=1):
    bpass = AIPSTask('bpass', version = aipsver)
    bpass.indata = uvdataset
    bpass.calsour[1] = srcname
    bpass.docalib = 1
    bpass.gainuse = clversion
    bpass.flagver = 1
    bpass.outver = 1
    bpass.solint = 0
    bpass.refant = refant
    bpass.bpassprm[5] = 1
    bpass.bpassprm[9] = 1
    bpass.bpassprm[10] = 1
    print "CL version is " + str(clversion) + ", calsour is " + str(bpass.calsour[1])
    bpass()

                                                                                                    
##### Selfcal using CALIB #################################
def run_calib(uvdataset, srcname, flux, snver, clver, refant=1):
    calib = AIPSTask('calib', version = aipsver)
    calib.indata = uvdataset
    calib.calsour[1] = srcname
    if flux is None:
        calib.smodel[1] = 1.0
        calib.smodel[2:] = [0]
    else:
        calib.smodel[1] = flux
        calib.smodel[2:] = [0]
    calib.docalib = 1
    calib.gainuse = clver
    calib.flagver = 1
    calib.doband = 1
    calib.bpver = 1
    calib.cmethod = 'DFT'
    calib.refant = refant
    calib.solint = 1
    calib.aparm[1] = 3
    calib.aparm[7] = 6
    calib.soltype = 'L1R'
    calib.solmode = 'A&P'
    calib.snver = snver
    calib()
    
# Flag the calibrator data, if desired
if args.flagfile is not None:
    flag_data(caldata, args.flagfile, 1)

# Run FRING
snversion       = 1  # Assume no existing SN table 
clversion       = 1
run_fring(caldata)

# Calibrate
applySN(caldata, snversion, clversion, refant)
snversion += 1
clversion += 1

# Bandpass solution
bpversion = 1
run_bpass(caldata, args.sourcename, clversion, refant)

# Run selfcal
run_calib(caldata, args.sourcename, args.flux, snversion, clversion, refant)

# Calibrate
applySN(caldata, snversion, clversion, refant)
snversion += 1
clversion += 1

