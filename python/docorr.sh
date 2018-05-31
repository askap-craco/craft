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

#outdir=/data/TETHYS_1/ban115/correlations/virgo/mode3/SB1633/card3.extrapause/
outdir=/data/TETHYS_1/ban115/correlations/virgo/mode1/SB1679/allcards/
mkdir -p $outdir
chmod a+w $outdir

if [[ ! -w $outdir ]] ; then
    echo "$outdir not writable"
    exit 1
fi

for f in $@ ; do
    dlname=`basename $f`
    echo "DL NAME IS $dlname"
    for b in beam40 beam41; do
	tsp craftcor.py --parset $fcm --calcfile $calcfile  -i 1024 -o $outdir/${dlname}_${b}.fits $dlname/c*/$b/*.vcraft
    done
done

