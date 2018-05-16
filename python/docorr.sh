#!/bin/bash

#fcm=$1
#calcfile=$2
#dlname=`basename $3`
#outdir=$4

if [[ ! -e $fcm ]] ; then
    echo No FCM $fcm
    exit 1
fi

if [[ ! -e $calcfile ]] ; then
    echo No calcfile $calcfile
    exit 1
fi

outdir=/data/TETHYS_1/ban115/correlations/virgo/mode0/SB1633/card3.extrapause/
mkdir -p $outdir
chmod a+w $outdir

if [[ ! -w $outdir ]] ; then
    echo "$outdir not writable"
    exit 1
fi

for f in $@ ; do
    dlname=`basename $f`
    echo "DL NAME IS $dlname"
    for b in beam00 beam01; do
	tsp craftcor.py --parset $fcm --calcfile $calcfile  -i 1024 -o $outdir/${dlname}_${b}.fits `find $dlname -name '*.vcraft' | grep $b | grep -v co34`
    done
done

