#!/bin/bash

candfile=fredda_galaxy18.cand

for f in `find -L . -name 'ak*.hdr'` ; do
    fdir=`dirname $f`
    pushd $fdir > /dev/null

    if [[ ! -s $candfile || ! -s $candfile.0 || ! -s $candfile.1 || ! -s $candfile.2 || ! -s $candfile.3 ]] ; then
	rm -f $candfile
	echo $fdir
    fi

    popd > /dev/null
done
