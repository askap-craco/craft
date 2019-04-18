#!/bin/bash
# ---------------------------------------------------------------------------------
# Shell script to reconstruct a localized ASKAP FRB to its full time resolution (~3 ns)
# This script consists of 2 parts:
#   PART 1) Create high resolution 1D array in frequency domain with delay corrections
#   PART 2) Remove ripples generated at PFB stage, coherently dedisperse, and ifft back to time series

# NOTE : i=1, n=16384 takes a long time, so it is recommended to run parallel for each antenna and then sum up the time series at the very last step
# (for each antenna: max. time ~ 9 hrs, max. memory ~ 30g)
# maybe consider overlap-save method in the future
# ---------------------------------------------------------------------------------


# PART 1
## Specify one antenna number, or dont pass argument if you want to add all antennas
an=$1

args=()
## set array size and offset
i=1
n=8 #16384 # $n * 54 * 336 is the total length of the output array
args+=("-i $i")
args+=("-n $n")
offset=1874193
args+=("--offset $offset")

if [ $# != 0 ]
then
args+=("--an ${an}")
fi

## Calibration solution directories - change
basedr=/fred/oz002/users/hcho/craft/
### geometric delays
calcfile=${basedr}Calibration/FRB180924-calcfiles/SB6635_frbpos.im
### clock delays
fcm=${basedr}Calibration/FRB180924-calcfiles/fcm_release2_ak25mod.txt
### hardware delays
hwfile=${basedr}Calibration/SB6635_b18_neweop.hwdelays
### gain, bandpass delays
mir=${basedr}Calibration/0407_allfreq/20180924212734_call_beam37_i1024_f9.uvaver
coeff=${basedr}Calibration/aips_bp_polyfit_coeff8.npy
args+=("--calcfile $calcfile")
args+=("--parset $fcm")
args+=("--hwfile $hwfile")
args+=("--mirsolutions $mir")
args+=("--aips_c $coeff") # optional, comment out if you want to use MIRIAD bandpass corrections instead of AIPS

## Part 1 output directory (frequency domain) - change
f_outfile=./test_f.npy #${basedr}python/test_output/test_f_an${an}.npy
args+=("-o $f_outfile")

## VCRAFT file directory - change
f_vcraft=${basedr}python/voltages/FRB180924/ak**/beam37/*.vcraft

python craftcor.py ${args[@]} --tab $f_vcraft


# PART 2
DM=362.193
t_outfile=./test_t.npy #${basedr}python/test_output/test_t_an${an}.npy
fftlen=$(( $n*64 ))

python freq2time.py -f $f_outfile -d $DM -o $t_outfile -l $fftlen -q
