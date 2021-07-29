#!/bin/bash

#SBATCH -A bsd
#SBATCH -p batch
#SBATCH -t 0-00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem-per-cpu=2G
#SBATCH -J pyname
#SBATCH -o outdir/out.%J
#SBATCH -e outdir/err.%J

module load gromacs

export GMX_MAXBACKUP=-1;
export OMP_NUM_THREADS=24;

echo "begin job.."
echo $PWD


# Inputs
msdout="msd_nptmain"; outdir="msdanalysis"
nchains=py_nchains

if ! test -f "py_trajfile"; then
    printf "py_trajfile not found"
    exit 1
elif ! test -f "py_tprfile"; then
    printf "py_tprfile not found"
    exit 1
elif ! test -f "py_conffile"; then
    printf "initconf not found"
    exit 1
elif ! test -f "chaininp.inp"; then
    echo "chaininp.inp not found"
    exit 1
fi
	

# Make Index files
if ! test -f "chindx.ndx"; then 
    srun gmx select -f py_conffile -s py_tprfile -sf chaininp.inp -on chindx.ndx
fi
wait

# Compute Rg of chains
printf "Computing Rg of chains"

if [ ! -d ${outdir} ]; then
    mkdir -p ${outdir}
fi
wait

for (( chcnt_i = 0; chcnt_i < nchains-1; chcnt_i++ ))
do
    printf "${chcnt_i}" | srun gmx msd -f py_trajfile -s py_tprfile -n chindx.ndx -o ${msdout}_${chcnt_i}.xvg &
done
wait

mv ${msdout}_*.xvg ${outdir}
cp chainlist.dat ${outdir}
cp chindx.ndx ${outdir}

printf "End of Rg calculations.."
