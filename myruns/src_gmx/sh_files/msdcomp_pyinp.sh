#!/bin/bash

#SBATCH -A bsd
#SBATCH -p batch
#SBATCH -t 0-03:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G
#SBATCH -J pyname
#SBATCH -o outdir/out.%J
#SBATCH -e outdir/err.%J

module load PE-gnu/3.0
module load cuda/10.1
module load gromacs/2020.6

export OMP_NUM_THREADS=24;

echo "begin job.."
echo $PWD


# Inputs
msdout="msd_nptmain"; allresultdir="all_results"
nchains=py_nchains
	
# Make Index files
if ! test -f "chaininp.inp"; then
    echo "chaininp.inp not found"
    exit 1
fi

if ! test -f "chindx.ndx"; then 
    srun gmx select -f py_conffile -s py_tprfile -sf chaininp.inp -on chindx.ndx
fi
wait

# Search for tpr/trr files
if ! test -f "enermin.tpr"  && ! test -f "traj_npt_main.trr"; then
    printf "tpr/trr files not found"
    exit 1
fi

if ! test -f "traj_npt_main_nojump_100ps.trr"; then
    printf "0" | srun gmx trjconv -s enermin.tpr -f traj_npt_main.trr -dt 100 -pbc nojump -o traj_npt_main_nojump_100ps.trr
wait
fi


# Compute MSD of chains
printf "Computing MSD of chains"

mkdir -p ${allresultdir}

for (( chcnt_i = 0; chcnt_i < nchains-1; chcnt_i++ ))
do
    printf "${chcnt_i}" | srun gmx msd -f traj_npt_main_nojump_100ps.trr -s py_tprfile -n chindx.ndx -o ${msdout}_${chcnt_i}.xvg &
done
wait

mv ${msdout}_*.xvg ${allresultdir}
cp chainlist.dat ${allresultdir}
cp chindx.ndx ${allresultdir}

printf "End of MSD calculations.."
