#!/bin/bash
#$ -S /bin/bash  # shell
#$ -cwd  # current working directory
#$ -j y  # merge stdout/stderr into one file
#$ -V  # make current environment variables available
#$ -o artemis.out  # stdout stream

## qsub -t 1-3 -pe chost 1 -l h_rt=1:00:00,h_vmem=2G $yorozuya/toolkit/transients/qsub_artemis_list.sh <outputdir> <pix.list> <-1> <dt> <inj.npz>

unset DISPLAY

echo job $SGE_TASK_ID
echo host $HOSTNAME

SCRIPTDIR=/csr/mwa/lufeng/toolkit/transients

rundir=$1  # output directory (must already exist)
pixlist=$2
j=$3  # single template run (unless it's -1)
dt=$4  # duration shift in sec
inj=$5  # injection npz file

## in case I want to loop over duration from within the script
## number of images - 2 for number of durations to iterate
#nitems=`wc -l im.list | awk '{print $1}'`
#nitems=$(($nitems - 2))
#    for j in `seq 0 $nitems`

i=0
while read -r pixarr
do

    if [ ${j} -ne -1 ]; then
	### single template run
	if [ ${inj} = "none" ]; then
	    ### no injection
	    python ${SCRIPTDIR}/artemis.py -i im.list -b bm.list -s $pixarr -o ${rundir}/${SGE_TASK_ID}r${i}.fits -l -c ${j},0 -d ${dt}
	else
	    ### yes injection
	    python ${SCRIPTDIR}/artemis.py -i im.list -b bm.list -s $pixarr -o ${rundir}/inj_${SGE_TASK_ID}r${i}.fits -l -c ${j},0 -j -f ${inj} -d ${dt}
	fi
    else
	### many template run
	if [ ${inj} = "none" ]; then
	    ### no injection
	    python ${SCRIPTDIR}/artemis.py -i im.list -b bm.list -s $pixarr -o ${rundir}/${SGE_TASK_ID}r${i}.fits -d ${dt}
	else
	    ### yes injection
	    python ${SCRIPTDIR}/artemis.py -i im.list -b bm.list -s $pixarr -o ${rundir}/inj_${SGE_TASK_ID}r${i}.fits -j -f ${inj} -d ${dt}
	fi
    fi

## in case I ever want to loop over start times instead of durations
#    mkdir output_t${j}
#    python ${SCRIPTDIR}/artemis.py -i im.list -b bm.list -s $pixarr -o output_t${j}/inj_${SGE_TASK_ID}${i}.fits -l -c 47,${j} -j -f single_injparam.npz

    i=$(($i + 1))

done < "${SGE_TASK_ID}_$pixlist"
