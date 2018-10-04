#!/bin/bash

w=1

infits=ak01-2_c4_f3_i16.fits
outfits=ak01-2_c4_f3_i16_w${w}.fits
rm -f $outfits

uvfrbwt.py $infits -o $outfits -w $w -c snoopy.log.finf_t51 --weight-data
uvfits2fil.py $outfits --weight
uvfits2fil.py $infits

plot_allbeams.py ak01-2_c4_f3_i16_w${w}-ak01-ak01-auto.fil -d 361

