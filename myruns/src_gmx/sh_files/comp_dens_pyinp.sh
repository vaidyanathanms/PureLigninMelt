#!/bin/bash

#SBATCH -A bsd
#SBATCH -p burst
#SBATCH -t 0-01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH -J py_jname
#SBATCH -o outdir/out.%J
#SBATCH -e outdir/err.%J

module load PE-gnu/3.0
module load gromacs/2020.6
module load vmd

export OMP_NUM_THREADS=32;

headdir=${PWD}
enerfile_1="ener_npt_main.edr"; tprfile_1="npt_main.tpr"
densfile_1="dens_npt.xvg"

allresultdir="all_dens"	
for temp_i in {py_Tinit..py_Tfin..py_dT}
    do
	tempdir="./T_${temp_i}"
	if [ ! -d "${tempdir}" ]; then
	    echo "${tempdir} not found"
	    continue
	else
	    echo "Analyzing density data in ${tempdir}.."
	    cd ${tempdir}
	fi

	if ! test -f "$enerfile_1" || ! test -f "$tprfile_1" ; then
	    echo "$enerfile_1 or $tprfile_1 not found"
	else
	    printf "22\n0" | srun gmx energy -f $enerfile_1 -s $tprfile_1 -o $densfile_1 &
	fi
	wait

        mkdir -p ${allresultdir}
        mv dens_*xvg ${allresultdir}
	wait

	cd ${headdir}
    done

