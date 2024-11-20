#!/bin/bash

THREADS=16

mkdir -p ./groups/
for peak in ./bed/PC3_WT_FOXA2_RECHIP_* ; do
    pname=$(echo $peak | awk -F '/' '{ gsub(".bed", "" ); print $NF}')
    echo $pname
    if [ -d ./groups/${pname} ] ; then
        echo "${pname} already processed "
    else
        findMotifsGenome.pl $peak hg38 ./groups/${pname}  -size 200 -p ${THREADS}
        findMotifsGenome.pl $peak hg38 ./groups/${pname}_given  -size given -p ${THREADS}

    fi
done
