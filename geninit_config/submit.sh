#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --job-name=geninp
#SBATCH -o ./geninp_%j.out
#SBATCH -e ./geninp_%e.err
#SBATCH -p testing
#SBATCH -A bsd

echo ${SLURM_SUBMIT_DIR}

module rm PE-gnu
module load PE-intel
module load python
module load namd
module load vmd

python initialize_dirs_for_runs.py inputsforpsfgen.inp
