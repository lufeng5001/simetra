#!/bin/bash

run=$1
hold=$2
j=$3

if [ $# -lt 3 ]; then
    echo "Usage: ./submit.sh <runname> <holdname> <suffix>"
    exit
fi

for (( i=1; i<=46; i++ ))
do

    /csr/mwa/lufeng/toolkit/transients/artbatch.sh ${run}_${i}${j}_config.txt $i ${hold}_${i}${j}_rr0

done
