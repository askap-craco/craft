#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=24

export ACES=/home/ban115/ACES
export CRAFT=/home/ban115/craft/craft/
export PATH=$CRAFT/cuda-fdmt/cudafdmt/src:$CRAFT/python:$PATH:/home/ban115/bin:$CRAFT/jobs/
export PYTHONPATH=$ACES/pythonlib:$CRAFT/python:$PYTHONPATH
export OMP_NUM_THREADS=24

module unload PrgEnv-cray
module unload gcc/4.9.0
module load python/2.7.10
module load astropy
module load pyephem



thedir=$1
thebase=`basename $thedir`
# Array = co or ak
array=co
pset=/home/ban115/fixheaders/$array/sbpars/
parset=$pset/${thebase}.parset
sblist=$pset/sblist.txt

echo "`date` Fixing SB in directory $thedir with base $thebase. Parset=$parset sblist=$sblist on host `hostname`"

if [[ ! -e $parset ]] ; then
    echo Parset $parset doesnt exist
    exit 1
fi

if [[ ! -e $sblist ]] ; then
    echo SBLIST $sblist doesnt exist
    exit 1
fi

find $thedir -name "${array}*.hdr" | xargs fix_headers.py --parset-dir $pset  --sblist  $sblist --fix
find $thedir -name "${array}*.hdr.v2" | xargs summarise_scans.py --write-meta-files --outfile /home/ban115/archive/v9/$thebase.json

echo "`date` Finished fixing SB in directory $thedir with base $thebase. Parset=$parset sblist=$sblist on host `hostname`"
