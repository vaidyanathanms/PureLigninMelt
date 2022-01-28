#!/bin/bash

#SBATCH -A bsd
#SBATCH -p burst
#SBATCH -t 0-02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH -J py_jobname
#SBATCH -o outdir/out.%J
#SBATCH -e outdir/err.%J

module load PE-gnu/3.0
module load gromacs/2020.6
module load vmd

export OMP_NUM_THREADS=32;

headdir=${PWD}
enerfile_1="ener_npt_berendsen.edr"; tprfile_1="npt_berendsen.tpr"
enerfile_2="ener_npt_main.edr"; tprfile_2="npt_main.tpr"
densfile_1="dens_berendsen.xvg"; densfile_2="dens_npt.xvg"

for run_i in {1..1..1} 
do
    rundir="."
    if [ ! -d "${rundir}" ]; then
	echo "${rundir} not found"
	continue
    fi
	
    for temp_i in {py_Tinit..py_Tfin..py_dT}
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
	    printf "22\n0" | srun gmx_mpi energy -f $enerfile_1 -s $tprfile_1 -o $densfile_1 &
	fi

	if ! test -f "$enerfile_2" || ! test -f "$tprfile_2" ; then
	    echo "$enerfile_1 or $tprfile_1 not found"
	else
	    printf "22\n0" | srun gmx_mpi energy -f $enerfile_2 -s $tprfile_2 -o $densfile_2 &
	fi
	wait

	cd ${headdir}
    done
done
