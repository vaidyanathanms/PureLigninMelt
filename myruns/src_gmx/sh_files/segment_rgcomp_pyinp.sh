#!/bin/bash

#SBATCH -A bsd
#SBATCH -p batch
#SBATCH -t 0-10:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=2G
#SBATCH -J py_jname
#SBATCH -o outdir/out.%J
#SBATCH -e outdir/err.%J

module load PE-gnu/3.0
module load cuda/10.1
module load gromacs/2020.6
module load vmd

export GMX_MAXBACKUP=-1;
export OMP_NUM_THREADS=32;

echo "begin job.."
echo $PWD


# Inputs
rgout="traj_npt_main.trr"; outdir="rganalysis"
nchains=py_nchains


if [ ! test -f "enermin.tpr"] && [ ! test -f "traj_npt_main_nojump_100ps.trr"]; then
    printf "tpr/trr files not found"
    exit 1
fi

#if ! test -f "traj_npt_main_nojump_100ps.trr"; then
#    printf "0" | srun gmx trjconv -s enermin.tpr -f traj_npt_main.trr -dt 100 -pbc nojump -o traj_npt_main_nojump_100ps.trr
#wait
#fi

if ! test -f "../init_files/L.psf"; then
    printf "psf file not found"
    exit 1
fi
	

# Compute segmental Rg of chains
printf "Computing segmental Rg of chains"

if [ ! -d ${outdir} ]; then
    mkdir -p ${outdir}
fi
wait

# Segmental Rg computation command using tcl script
vmd -dispdev text -e calc_seg_rg.tcl

mv rg_allsegs.dat ${outdir}
mv RgvsN.dat ${outdir}

printf "End of segmental Rg calculations.."



