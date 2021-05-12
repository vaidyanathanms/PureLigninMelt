#!/bin/bash

headdir=${PWD}
tprfile_old="npt_main.tpr"; enerfile_old="ener_npt_main.edr"; trajfile_old="traj_npt_main.trr"
tprfile_new="npt_main.tpr" # dont change the name
extendtime=100000 #ps

for run_i in {1..4..1}
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

	if ! test -f "$tprfile_old"; then
	    echo "$tprfile_old/$enerfile_old/$trajfile_old not found in $dirpath"
	else
	    jsrun -X 1 -n 1 -c 1 -a 1 -g 1 --launch_distribution plane:1 -b packed:1 gmx_mpi convert-tpr -s $tprfile_old -o $tprfile_new -extend ${extendtime}
	fi

	cd ${headdir}
    done
done
