#!/bin/bash

echo Doing $0 $@ scan=SB0${SBID}/${SCANNAME}/ICS/${CID} alias="${SB_ALIAS}"
fname=$1
fredfof_ics.py $@
plot_fredda_cand.py ${fname}.fof -o ${fname}.fof.png
#postslack.py -t 'Scan summary' scan=SB0${SBID}/${SCANNAME}/ICS/${CID} alias=""${SB_ALIAS}"" -i ${fname}.fof.png



postslack.py -t 'Scan summary' scan=SB${SBID}/${SCANNAME}/ICS/${CAPTURE_NAME} alias="${SB_ALIAS}" -i ${fname}.fof.png
#postslack.py -t 'Scan Complete'
if [ ! -f ${fname}.fof ]; then
       postslack.py -t 'Completed scan details' scan=SB${SBID}/${SCANNAME}/ICS/${CAPTURE_NAME} alias="${SB_ALIAS}"
fi


# select a few candidates for plotting

icsdir=/data/TETHYS_1/craftop/data/SB${SBID}/${SCANNAME}/ICS/${CAPTURE_NAME}/
echo  $icsdir


echo  "fof name ${fname}.fof"
ls ${fname}.fof
wc ${fname}.fof


sort -gk 1 ${fname}.fof  | awk '{if ( ($1 > 9.5) && ($3 > 10) && ($4 < 10) && ($6 >100)) print $1,$2,$3,$4,$5,$6,$7,$8,0}' | awk '{if (NR < 5) print $0}' > ${fname}_good.list

echo "size of good file"
wc ${fname}_good.list
ngood=`wc -l ${fname}_good.list | awk '{print $1}'`


echo "ngood is $ngood"
 
for igood in `seq 1 $ngood`;  do
    echo "# S/N, sampno, secs from file start, boxcar, idt, dm, beamno, mjd, latency_ms" > ${fname}_${igood}.cand
    btmp=`sed "${igood}q;d" ${fname}_good.list | awk '{print $7}'`
    if [ $btmp -lt 10 ]; then
        beamno=0${btmp}
    else
        beamno=$btmp;
    fi    
    echo $btmp

    sed "${igood}q;d" ${fname}_good.list >> ${fname}_${igood}.cand
    #cat ${fname}_${igood}.cand | /data/TETHYS_1/bha03n/test/scripts/postsnoo.py     --sbid $SBID --scan $SCANNAME --alias "$SB_ALIAS" --ant "$ANTENNAS_COMMA" --cid $CAPTUREID --text "Offline (fredda fof)  ICS FRB  $FOOTPRINT_NAME $FIELD_NAME" 
    echo ${icsdir}
    /data/TETHYS_1/bha03n/test/scripts/filplt_zero_cand_tscrunch.py -t ${INT_TIME} -c ${fname}_${igood}.cand   ${icsdir}/*.fil -f $fname
    #postslack.py -i /data/TETHYS_1/craftop/data/SB${SBID}/${SCANNAME}/ICS/${CAPTURE_NAME}/*${beamno}.png -t "candidate plot"
        
done


    
#postslack.py -t 'Scan Complete' 
