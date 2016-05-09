#!/bin/bash
#$ -N ioconversion  # job name
#$ -S /bin/bash  # shell
#$ -cwd  # current working directory
#$ -j y  # merge stdout/stderr into one file
#$ -V  # make current environment variables available
#$ -o ioconversion.out  # stdout stream

## qsub -t 11 -pe chost 1 -l h_rt=2:00:00,h_vmem=3G $yorozuya/toolkit/transients/qsub_ioconversion.sh pix.list

unset DISPLAY

SCRIPTDIR=/csr/mwa/lufeng/toolkit/transients

pixlist=$1

while read -r pixarr
do

    python ${SCRIPTDIR}/artemis.py -i im.list -b bm.list -s $pixarr -w

done < "${SGE_TASK_ID}_$pixlist"
