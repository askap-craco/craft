#!/bin/bash

# Links the directory hierarchy 
root=/scratch2/askap/askapops/craft/

for f in $@ ; do
    f=`dirname $f`
    mainpath=$root/$f
    if [[ ! -d $mainpath ]] ; then
	echo $mainpath does not exist
	exit 1
    fi

    if [[ ! -d $f ]] ; then
	mkdir $f
    fi

    # make directories - this makes me hate myself 
    find $mainpath -type d | sed "s!$root!!" | sed "s!/!!" | xargs mkdir -p
    #find $mainpath -type d

    for d in `find $mainpath -type d -not -empty` ; do
	relpath=`echo $d | sed "s!$root!!"`
	echo "Relpath is $relpath"
	pushd `dirname $relpath`
	echo ln -s $d/* .
	popd

    done

done
