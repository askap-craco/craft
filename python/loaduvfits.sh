#!/bin/bash

for f in $@ ; do
    f=`echo $f | sed s/.fits//`

    if [[ ! -e $f.mir ]] ; then
	fits in=$f.fits out=$f.mir op=uvin
	prthd in=$f.mir
	uvflag vis=$f.mir flagval=unflag select=-auto
	uvspec vis=$f.mir select=ant\(1\) interval=1 axis=freq,real log=$f.uvspec.real
	uvspec vis=$f.mir select=ant\(1\) interval=1 axis=freq,imag log=$f.uvspec.imag
	mfcal vis=$f.mir interval=1
	gpplt vis=$f.mir log=$f.gains.real yaxis=real
	gpplt vis=$f.mir log=$f.gains.imag yaxis=imag
	gpplt vis=$f.mir log=$f.bandpass.real yaxis=real options=bandpass
	gpplt vis=$f.mir log=$f.bandpass.imag yaxis=imag options=bandpass
	
    fi

    if [[ ! -e $f.uvaver ]] ; then
	uvaver vis=$f.mir out=$f.uvaver line=chan,336,1,6 options=nopass,nocal,npol
	mfcal vis=$f.uvaver interval=1
	uvspec vis=$f.uvaver axis=freq,phase device=$f.uvaver.uvspec.png/png interval=1 nxy=3,4

	uvplt vis=$f.uvaver axis=time,phase device=$f.uvaver.uvplt.png/png nxy=3,4
    fi
	
done

d=`dirname $0`

