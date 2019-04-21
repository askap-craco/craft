#!/bin/bash

for f in $@ ; do
    f=`echo $f | sed s/.fits//`

    if [[ ! -e $f.mir ]] ; then
	fits in=$f.fits out=$f.mir op=uvin
	prthd in=$f.mir
	uvflag vis=$f.mir flagval=unflag select=-auto
	uvspec vis=$f.mir select=ant\(1\) interval=1 axis=freq,real log=$f.uvspec.real
	uvspec vis=$f.mir select=ant\(1\) interval=1 axis=freq,imag log=$f.uvspec.imag
	calibrate.sh $f.mir
	
    fi

    if [[ ! -e $f.uvaver ]] ; then
	#uvaver vis=$f.mir out=$f.uvaver line=chan,336,1,6 options=nopass,nocal,nopol interval=1
	uvaver vis=$f.mir out=$f.uvaver line=chan,48,1,42 options=nopass,nocal,nopol interval=1
	uvflag vis=$f.uvaver flagval=unflag select=-auto
	uvflag vis=$f.uvaver flagval=flag select=amp\(1.2\)
	#uvflag vis=$f.uvaver flagval=flag select=uvrange\(4,1e6\) line=chan,100,200
	uvlin vis=$f.uvaver out=$f.uvlin order=9 mode=cont options=nopass,nopol,nocal
	calibrate.sh $f.uvaver
	calibrate.sh $f.uvlin
	uvspec vis=$f.uvaver axis=freq,phase device=$f.uvaver.uvspec.png/png interval=1 nxy=3,4

	uvplt vis=$f.uvaver axis=time,phase device=$f.uvaver.uvplt.png/png nxy=3,4
    fi
	
done

d=`dirname $0`

