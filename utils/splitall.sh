#!/bin/bash

for f in `find . -name '*.fil'` ; do 
    dirname=`dirname $f`
    fname=`basename $f`
    pushd $dirname
    splitfil.py $fname
    popd
done
