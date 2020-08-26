#!/bin/bash

function dosim() {
    ( set -o posix ; set) > /tmp/myvars.before
    frb_idm=$1
    t0=$2 # start time in samples
    frb_amp=$3
    frb_sn=$4
    frb_relpos=$5

    # telescope and processing parameters
    ncu=2
    fch1=0.736 # GHz
    nt=16
    nd=4
    nchan=256 # number of channels
    tint=1.728
    desired_amp=500 # desired amplitude a the output of the FFT
    threshold=400
    antfile=/data/craco/ban115/craft/python/askap-ak1-ak20.ant

    # Calculated parameters
    frb_relpos_name=$( echo "$frb_relpos" | sed s/,//)
    tstart=$(echo "scale=6; $t0 * $tint " | bc) # start time in milliseconds
    nant=$(wc -l $antfile | cut -f 1 -d ' ')
    nbl=$( echo "scale=6; $nant * ($nant -1) / 2" | bc) # number of baselines
    name=frb_d${frb_idm}_lm${frb_relpos_name}_t${t0}_sn_${frb_sn}_nt${nt}_nant${nant}
    fits=${name}.fits

    # Expected output of the FFT amplitude
    # Not sure where the factor of 2 is coming from
    raw_image_amp=$( echo "scale=6; $nbl * $nchan * $frb_amp * 2" | bc )
    scale=$( echo "scale=6; $desired_amp / $raw_image_amp " | bc )

    echo $ncu $nt $nd $frb_amp $nant $nbl $nchan $desired_amp $raw_image_amp scale=$scale
    mkdir $name
    pushd $name
    cp $antfile .
    ( set -o posix ; set) > /tmp/myvars.after
    diff /tmp/myvars.before /tmp/myvars.after > varchanges.txt
    cp /tmp/myvars.after allvars.txt

    cmd="uvfrbsim.py --fch1 $fch1 --nchan $nchan --antfile $antfile --tint $tint --duration $nt --frb_idm $frb_idm --frb_amp $frb_amp --frb_sn $frb_sn --frb_relpos $frb_relpos -o $fits"
    echo "Running $cmd"
    $cmd

    craco_fdmt_krnl.py --nt $nt --ndm $nd --format raw --nfftcu $ncu --output-scale $scale $fits
    craco_img_krnl.py --ndm $nd --uvgrid $fits.uvgrid.txt --nfftcu $ncu --nt $nt  $fits.ndm${nd}_nt${nt}.b0.uvdata.raw --threshold $threshold

    popd
}

#    frb_idm=$1
#    t0=$2 # start time in samples
#    frb_amp=$3
#    frb_sn=$4
#    frb_relpos=$5

dosim 0 0 1 inf 0,0 # Canonical:  unit FRB at center at 0 DM
dosim 3 0 1 inf 0,0 # canonical + idm=3
dosim 0 4 1 inf 0,0 # canonical + toff=4 samlpes
dosim 0 0 2 inf 0,0 # canonical + amplitude=2
dosim 0 0 1 inf 100,200 # canonical + position offset = 100" in ra and 200" in dec
dosim 0 0 1 10 0,0 # Canonical + S/N = 10
dosim 3 4 2 inf 300,400 # Everything with infinite S/N
dosim 3 4 2 10 300,400 # Everything with S/N=10
