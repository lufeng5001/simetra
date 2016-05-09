#!/bin/bash

## This script submits artemis jobs in batch with a unique run name and keeps track of rerun numbers

# user input
program=`basename $0`
config=$1
pixarray=$2
holdname=$3

# usage
if [ $# -lt 2 ]; then
    echo "Usage:" ${program} "<config.txt> <pixarray> <holdname>"
    exit
fi

# read artemis config file
source ${config}

# create run directory and keep track of rerun number
if [ ! -d $RUNNAME ]; then
    mkdir -v $RUNNAME  # first run
    rrnum=0
else
    read -r rrnum < "$RUNNAME/rerun"  # rerun
    rrnum=$(($rrnum + 1))
fi
echo $rrnum > $RUNNAME/rerun  # log current run number for future reference

# submit jobs
if [[ $pixarray == *-* ]]; then
    IFS='-' read -ra pixarr <<< "${pixarray}"  # parse num1-num2 into an array; create array range between num1 and num2
    for i in $(seq ${pixarr[0]} ${pixarr[1]}); do
	if [ -z ${holdname} ]; then
	    qsub -N ${RUNNAME}_rr${rrnum} -t $i -pe chost 1 -l h_rt=48:00:00,h_vmem=3G $yorozuya/toolkit/transients/qsub_artemis_list.sh $RUNNAME $PIXFILE $TEMPLATE $DURATION $INJFILE
	else
	    qsub -hold_jid ${holdname} -N ${RUNNAME}_rr${rrnum} -t $i -pe chost 1 -l h_rt=48:00:00,h_vmem=3G $yorozuya/toolkit/transients/qsub_artemis_list.sh $RUNNAME $PIXFILE $TEMPLATE $DURATION $INJFILE
	fi
    done
    unset IFS
else
    IFS=','
    for i in ${pixarray}; do  # split num1,num2,...,num3 into an array
	qsub -N ${RUNNAME}_rr${rrnum} -t $i -pe chost 1 -l h_rt=48:00:00,h_vmem=3G $yorozuya/toolkit/transients/qsub_artemis_list.sh $RUNNAME $PIXFILE $TEMPLATE $DURATION $INJFILE
    done
    unset IFS
fi
