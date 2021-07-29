#!/bin/bash

headdir=${PWD}
initfile="initconf.gro"
tprfile="npt_main.tpr"
trajfile="traj_npt_main.trr"
rginp="rginp.inp"
rgout="rg_nptmain"
nchains=22

for run_i in {6..6..1}
do
    rundir="./run_${run_i}"
    if [ ! -d "${rundir}" ]; then
	echo "${rundir} not found"
	continue
    fi
	
    for temp_i in {400..400..20}
    do
	tempdir="/T_${temp_i}"
	dirpath="${rundir}${tempdir}"
	outdir="rgout_dir"
	if [ ! -d "${dirpath}" ]; then
	    echo "${dirpath} not found"
	    continue
	else
	    echo "Analyzing Rg data in ${dirpath}.."
	    cd ${dirpath}
	    echo $PWD
	fi

	if ! test -f "$tprfile" || ! test -f "$trajfile" || ! test -f "$rginp" || ! test -f $initfile; then
	    echo "$trajfile or $tprfile or $rginp not found"
	    cd ${headdir}
	    continue
	fi
	
	if [ ! -d "${outdir}" ]; then
 	     mkdir ${outdir}
	fi

	echo ${dirpath}

	# Make Index files

	jsrun -X 1 -n 1 -c 1 -a 1 -g 1 --launch_distribution plane:1 -b packed:1 gmx_mpi select -f $initfile -s $tprfile -sf $rginp -on rgindx.ndx	
	wait
	
	for chcnt_i in {0..21..1}
	do
	    echo ${chcnt_i}
	    rgfyle="$rgout_${chcnt_i}.xvg"
	    printf "${chcnt_i}" | jsrun -X 1 -n 1 -c 1 -a 1 -g 1 --launch_distribution plane:1 -b packed:1 gmx_mpi gyrate -f $trajfile -s $tprfile -n rgindx.ndx -o $rgfyle &
	done
	wait

	mv "$rgout_*.xvg" ${outdir}
	cd ${headdir}

    done
done
