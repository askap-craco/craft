#!/bin/bash -f

rm solved_fpgadelays.txt

for fpga in "$@"
do
    pushd $fpga
    echo "$fpga"
    difx2fits craftfrbD2D
    fpgadelay.py CRAFTFR.0.bin0000.source0000.FITS -o ../solved_fpgadelays.txt
    popd
done
