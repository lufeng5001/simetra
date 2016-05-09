#!/bin/bash

# This script generates injection files for each pix.list drawn from the uniform distribution.

for (( i=1; i<=46; i++ ))
do

    for j in a b c
    do
	python $yorozuya/toolkit/transients/injection.py -i /nfs/mwa-09/r1/128T_pipeline2/2013_09_02/1062173304/clean_fullband_1062173304_XX.fits,/nfs/mwa-13/d1/128T_pipeline5/2013_11_29/1069764984/clean_fullband_1069764984_XX.fits -o inj_${i}${j}.npz -n 86 -a uniform -d uniform -p uniform -m 0.15
	#python $yorozuya/toolkit/transients/injection.py -i /nfs/mwa-09/r1/128T_pipeline2/2013_09_02/1062173304/clean_fullband_1062173304_XX.fits,/nfs/mwa-13/d1/128T_pipeline5/2013_11_29/1069764984/clean_fullband_1069764984_XX.fits -o inj_fred_${i}${j}.npz -n 86 -a uniform -d uniform -p uniform -m 1 -t fred
    done
done
