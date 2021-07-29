#!/bin/bash

headdir=${PWD}
enerfile_1="ener_npt_berendsen.edr"; tprfile_1="npt_berendsen.tpr"
enerfile_2="ener_npt_main.edr"; tprfile_2="npt_main.tpr"
densfile_1="dens_berendsen.xvg"; densfile_2="dens_npt.xvg"

for run_i in {6..6..1}
do
    rundir="./run_${run_i}"
    if [ ! -d "${rundir}" ]; then
	echo "${rundir} not found"
	continue
    fi
	
    for temp_i in {300..500..20}
    do
	tempdir="/T_${temp_i}"
	dirpath="${rundir}${tempdir}"
	if [ ! -d "${dirpath}" ]; then
	    echo "${dirpath} not found"
	    continue
	else
	    echo "Analyzing density data in ${dirpath}.."
	    cd ${dirpath}
	fi

	if ! test -f "$enerfile_1"  || ! test -f "$tprfile_1" ; then
	    echo "$enerfile_1 or $tprfile_1 not found"
	else
	    printf "22\n0" | jsrun -X 1 -n 1 -c 1 -a 1 -g 1 --launch_distribution plane:1 -b packed:1 gmx_mpi energy -f $enerfile_1 -s $tprfile_1 -o $densfile_1 &
	fi

	if ! test -f "$enerfile_2" || ! test -f "$tprfile_2" ; then
	    echo "$enerfile_1 or $tprfile_1 not found"
	else
	    printf "22\n0" | jsrun -X 1 -n 1 -c 1 -a 1 -g 1 --launch_distribution plane:1 -b packed:1 gmx_mpi energy -f $enerfile_2 -s $tprfile_2 -o $densfile_2 &
	fi
	wait

	cd ${headdir}
    done
done
