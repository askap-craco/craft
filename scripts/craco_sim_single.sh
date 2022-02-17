#!/bin/bash

function param {
    name=$1
    shift
    echo $name $@ >> params.txt
    echo export $name=$@ >> params.env
}

if [[ $# != 5 ]] ; then
    echo "Need 5 arguments"
    echo "Usage: $0 DM T0 AMP SN RELPOS"
    echo "  DM - DM in PC/CM3"
    echo "  T0 - startign time in samples"
    echo "  AMP - Amplitude of output FRB"
    echo "  SN - Desired Signal to noise"
    echo "  RELPOS - comma separated ra_arcsec,dec_arcsec"
    exit 1
fi

( set -o posix ; set) > /tmp/myvars.before
frb_dm=$1
t0=$2 # start time in samples
frb_amp=$3
frb_sn=$4
frb_relpos=$5

# telescope and processing parameters
ncu=2
fch1=0.736 # GHz
nt=16
nd=8
nchan=256 # number of channels
tint=1.728
desired_amp=500 # desired amplitude a the output of the FFT
threshold=10
ncin=32
ndout=8
#antfile=`dirname $0`/askap-ak1-ak30.ant
#antfile=`dirname $0`/askap-ak1-ak3.ant
antfile="askap-ak1-ak24.ant"
antfile=$(antenna_locations $antfile)
echo Antfile is $antfile


# Calculated parameters
frb_relpos_frbname=$( echo "$frb_relpos" | sed s/,//)
frb_tstart=$(echo "scale=6; $t0 * $tint " | bc) # start time in milliseconds
nant=$(wc -l $antfile | awk '{print $1}')
nbl=$( echo "scale=6; $nant * ($nant -1) / 2" | bc) # number of baselines
frbname=frb_d${frb_dm}_t${t0}_a${frb_amp}_sn${frb_sn}_lm${frb_relpos_frbname}
fits=${frbname}.fits
sim_method=mkfrb_fdmt

# Expected output of the FFT amplitude
# Not sure where the factor of 2 is coming from
raw_image_amp=$( echo "scale=6; $nbl * $nchan * 2" | bc )
scale=$( echo "scale=6; $desired_amp / $raw_image_amp " | bc )

echo ncu=$ncu nt=$nt nd=$nd frb_amp=$frb_amp nant=$nant nbl=$nbl nchan=$nchan desired_amp=$desired_amp raw_amp=$raw_image_amp scale=$scale
mkdir $frbname
cp $antfile $frbname
pushd $frbname

antfile2=`basename $antfile`

( set -o posix ; set) > /tmp/myvars.after
diff /tmp/myvars.before /tmp/myvars.after > varchanges.txt
cp /tmp/myvars.after allvars.txt

rm -f params.txt
rm -f params.env

param rundate `date` 
param frb_dm $frb_dm 
param frb_t0 $t0 
param frb_amp $frb_amp
param frb_sn $frb_sn 
param frb_relpos $frb_relpos 
param frb_tstart $frb_tstart 
param ncu $ncu 
param fch1 $fch1 
param nt $nt
param nd $nd 
param nchan $nchan 
param tint $tint 
param desired_amp $desired_amp 
param threshold $threshold 
param antfile $antfile 
param frb_replos_frbname $frb_relpos_frbname 
param nant $nant 
param nbl $nbl 
param frbname $frbname 
param fits $fits 
param raw_image_amp $raw_image_amp 
param scale $scale
param ncin $ncin
param sim_method $sim_method

cmd="uvfrbsim --fch1 $fch1 --nchan $nchan --antfile $antfile2 --tint $tint --duration $nt --frb_idm $frb_dm --frb_amp $frb_amp --frb_sn $frb_sn --frb_relpos $frb_relpos --frb_tstart $frb_tstart -o $fits --sim-method $sim_method"
echo "Running $cmd"
#$cmd --show
$cmd

#craco_fdmt_krnl.py --nt $nt --ndm $nd --format raw --nfftcu $ncu --output-scale $scale $fits
#craco_img_krnl.py --ndm $nd --uvgrid $fits.uvgrid.txt --nfftcu $ncu --nt $nt  $fits.ndm${nd}_nt${nt}.b0.uvdata.raw --threshold $threshold

cmd="craco_pipeline --nt $nt --ndm $nd --ncin $ncin --ndout $ndout --uv $fits --boxcar-weight sqrt --threshold $threshold --save"

echo "Running $cmd"
$cmd

popd

