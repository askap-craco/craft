#!/bin/bash

#fcm=$1
#calcfile=$2
#dlname=`basename $3`
#outdir=$4

#suffix=$1
#shift

if [[ ! -e $fcm ]] ; then
    echo No FCM $fcm
    exit 1
fi

if [[ ! -e $calcfile ]] ; then
    echo No calcfile $calcfile
    exit 1
fi

#outdir=/data/TETHYS_1/ban115/correlations/virgo/mode3/SB1633/card3.extrapause/
#outdir=/data/TETHYS_1/ban115/correlations/virgo/mode3/SB1679/allcards/nfft4/

outdir=$corrdir/$suffix

echo "WRITING DATA TO $outdir"

mkdir -p $outdir
chmod a+w $outdir
ln -s $outdir/voltages `pwd`

if [[ ! -w $outdir ]] ; then
    echo "$outdir not writable"
    exit 1
fi

ismall=16

for f in $@ ; do
    dlname=`basename $f`
    echo "DL NAME IS $dlname"
    beams=`ls -d $f/*/beam37 |  awk -F / '{print $3}' | sort | uniq`
    echo "BEAMS ARE $beams"
    for b in $beams; do
	b=`basename $b`
	echo "beam is $b"
	tsp craftcor.py --parset $fcm --calcfile $calcfile  -i ${ismall} -o $outdir/${dlname}_c4_f3_${b}_i${ismall}.fits $dlname/ak??/$b/*c4_f3.vcraft --fft-size=1

	for c in {1..7} ; do
	    tsp craftcor.py --parset $fcm --calcfile $calcfile  -i 1024 -o $outdir/${dlname}_c${c}_${b}.fits $dlname/ak??/$b/*c${c}*.vcraft --fft-size=1
	done
    done
done

