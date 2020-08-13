#!/bin/bash


fits=frb_d0_lm0_nt16_nant24.fits
ncu=2
nt=16
nd=1
frbamp=1 # amplitude in uvfrbsim
nant=24
nbl=$( echo "scale=6; $nant * ($nant -1) / 2" | bc) # number of baselines
nchan=256 # number of channels
desired_amp=12 # desired amplitude a thte output of the FFT

# Expected output of the FFT amplitude
# Not sure where the factor of 2 is coming from
raw_image_amp=$( echo "scale=6; $nbl * $nchan * $frbamp * 2" | bc )
scale=$( echo "scale=6; $desired_amp / $raw_image_amp " | bc )

echo $ncu $nt $nd $frbamp $nant $nbl $nchan $desired_amp $raw_image_amp $scale

./uvfrbsim.py --antfile askap-ak1-ak24.ant --duration $nt -o $fits
./craco_fdmt_krnl.py --nt $nt --ndm $nd --format raw --nfftcu $ncu --output-scale $scale $fits
./craco_img_krnl.py --ndm $nd --uvgrid $fits.uvgrid.txt --nfftcu $ncu --nt $nt  $fits.ndm${nd}_nt${nt}.b0.uvdata.raw
