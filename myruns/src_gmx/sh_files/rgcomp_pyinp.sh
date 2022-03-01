#!/bin/bash

#SBATCH -A bsd
#SBATCH -p batch
#SBATCH -t 0-01:00:00
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

export GMX_MAXBACKUP=-1;
export OMP_NUM_THREADS=32;

echo "begin job.."
echo $PWD


# Inputs
rgout="rg_nptmain"; allresultdir="all_results"
eigout="eig_nptmain"
nchains=py_nchains

# Make Index files
if ! test -f "rgchaininp.inp"; then
    echo "rgchaininp.inp not found"
    exit 1
fi

if ! test -f "chindx.ndx"; then 
    srun gmx select -f py_conffile -s py_tprfile -sf rgchaininp.inp -on chindx.ndx
fi
wait

# Search for 100 ps interval file. If not present, create it
if ! test -f "traj_npt_main_nojump_100ps.trr"; then
    if ! test -f "npt_main.tpr"; then
	printf "tpr/trr files not found"
	exit 1
    else
	printf "0" | srun gmx trjconv -s npt_main.tpr -f traj_npt_main.trr -dt 100 -pbc nojump -o traj_npt_main_nojump_100ps.trr
	wait
fi

# Compute Rg/eigenvalues of chains
printf "Computing Rg of chains"

mkdir -p ${allresultdir}
wait

for (( chcnt_i = 0; chcnt_i <= nchains-1; chcnt_i++ ))
do
    printf "${chcnt_i}" | srun gmx gyrate -f traj_npt_main_nojump_100ps.trr -s py_tprfile -n chindx.ndx -o ${rgout}_${chcnt_i}.xvg &
    wait
    printf "${chcnt_i}" | srun gmx polystat -f traj_npt_main_nojump_100ps.trr -s py_tprfile -n chindx.ndx -nomw -pc -o ${eigout}_${chcnt_i}.xvg &
done
wait

mv ${rgout}_*.xvg ${allresultdir}
mv ${eigout}_*.xvg ${allresultdir}
cp chainlist.dat ${allresultdir}
cp chindx.ndx ${allresultdir}
cp rgchaininp.inp ${allresultdir}
printf "End of Rg calculations.."



