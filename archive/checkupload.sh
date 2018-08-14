#!/bin/bash

ldir=/group/askap/ban115/craft/archive-logfiles/
sbdir=/scratch2/askap/askapops/craft/co/
pdir=askap-craft/co
finfile=$ldir/archived_sbs.txt
cmpfiles=/group/askap/ban115/craft/archive/compares/

dodelete=0
while getopts "dh" opt  ; do
    case $opt in
	d)
	    dodelete=1
	    ;;
	h)
	    echo "$0 [-d] SB [SB [SB ...]]"
	    exit 0
	    ;;
	\?)
	    echo "Invalid option $OPTARG"
	    exit 1
	    ;;
    esac
done


cd $sbdir
shift $((OPTIND-1))

for f in $@ ; do
    f=`basename $f`
    echo "Checking $f"
    cmpfile="$cmpfiles/$f.compare"
    if [[ ! -d $f ]] ; then
	echo $f is not an SB
	continue
    fi

    pshellresult=`pshell "cd $pdir && compare $f" | grep -v meta | tee $cmpfile`
    pshellerr=$?
    echo "PSHELL returned code $pshellerr"
    if [[ $pshellerr != 0 ]] ; then
	echo "PSHELL FAILED. Quitting"
	exit 1
    fi

    nmissing=`awk '/=== Missing remote/{flag=1;next}/=== Compare complete/{flag=0}flag' < $cmpfile | wc -l`
    echo "Missing $nmissing" files
    if [[ $nmissing == 0 ]] ; then
	echo "$f No missing files"
    else
	echo "*** $f is missing $nmissing files"
    fi
    
    isfinished=`grep -x $f $finfile`
    echo "FINISHED file contains $isfinished"
    if [[ $isfinished == $f  ]] ; then
	echo "$f seems to have been successfully uploaded"
    else
	echo "$f *** IS NOT IN THE FINISHED FILE - PROBABLY NOT UPLOADED"
    fi

    if [[ $nmissing == 0 && $isfinished == $f ]] ; then
	echo "$f can be safely deleted"
	if [[ $dodelete == 1 ]] ; then
	    echo "Deleting filterbanks in $f"
	    find $f -name '*.fil' -delete
	    delret=$?
	    if [[ $delret == 0 ]] ; then
		echo "Filterbanks successfully deleted"
	    else
		echo "Delete failed with error code $delret"
		exit $delret
	    fi
	fi
    else
	echo "$f CANNOT be safely deleted"
    fi
done
