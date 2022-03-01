#!/bin/bash

#SBATCH -A bsd
#SBATCH -p batch
#SBATCH -t 0-05:00:00
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

export GMX_MAXBACKUP=-1;
export OMP_NUM_THREADS=24;

echo "begin job.."
echo $PWD

# Inputs
hbout="hb_nptmain"; hbout_dir="all_hbs"
nchains=py_nchains
printf "Computing hydrogen bonding"
mkdir -p ${hbout_dir}

# Compute inter/intra HB
for (( chcnt_i = 0; chcnt_i < nchains-1; chcnt_i++ ))
do
    srun gmx hbond -f py_trajfile -s py_tprfile -sf hb${chcnt_i}_chain.inp -o ${hbout}_${chcnt_i}.xvg
wait
done

mv ${hbout}_*.xvg ${hbout_dir}
mv hb*_chain.inp ${hbout_dir}
wait

# Compute FA-FA/FA-PCA/FA-G/FA-S
srun gmx hbond -f py_trajfile -s py_tprfile -sf hb_f.inp -o hb_f.xvg
mv hb_f.xvg ${hbout_dir}
mv hb_f.inp ${hbout_dir}
wait

printf "End of HB calculations.."
