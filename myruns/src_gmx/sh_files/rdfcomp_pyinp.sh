#!/bin/bash

#SBATCH -A bsd
#SBATCH -p burst
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
rdfout="rdf_nptmain"; rdfout_dir="all_rdfs"
nchains=py_nchains
printf "Computing RDF of chains"
mkdir -p ${rdfout_dir}

# Search for 100 ps interval file. If not present, create it
if ! test -f "traj_npt_main_nojump_100ps.trr"; then
    if ! test -f "npt_main.tpr"; then
	printf "tpr/trr files not found"
	exit 1
    else
	printf "0" | srun gmx trjconv -s npt_main.tpr -f traj_npt_main.trr -dt 100 -pbc nojump -o traj_npt_main_nojump_100ps.trr
	wait
    fi
fi

# Compute inter/intra RDF
for (( chcnt_i = 0; chcnt_i < nchains-1; chcnt_i++ ))
do
    srun gmx rdf -f traj_npt_main_nojump_100ps.trr -s py_tprfile -sf rdf${chcnt_i}.inp -o ${rdfout}_${chcnt_i}.xvg
wait
done

mv ${rdfout}_*.xvg ${rdfout_dir}
mv rdf*.inp ${rdfout_dir}
wait

# Compute H-H/H-G/H-S
srun gmx rdf -f traj_npt_main_nojump_100ps.trr -s py_tprfile -sf Hall.inp -o Hall.xvg
mv Hall.xvg ${rdfout_dir}
mv Hall.inp ${rdfout_dir}
wait

# Compute G-H/G-G/G-S
srun gmx rdf -f traj_npt_main_nojump_100ps.trr -s py_tprfile -sf Gall.inp -o Gall.xvg
mv Gall.xvg ${rdfout_dir}
mv Gall.inp ${rdfout_dir}
wait

# Compute S-H/S-G/S-S
srun gmx rdf -f traj_npt_main_nojump_100ps.trr -s py_tprfile -sf Sall.inp -o Sall.xvg
mv Sall.xvg ${rdfout_dir}
mv Sall.inp ${rdfout_dir}
wait

# Compute FA-FA/FA-PCA/FA-G/FA-S
srun gmx rdf -f traj_npt_main_nojump_100ps.trr -s py_tprfile -sf Fall.inp -o Fall.xvg
mv Fall.xvg ${rdfout_dir}
mv Fall.inp ${rdfout_dir}
wait

# Compute PCA-FA/PCA-PCA/PCA-G/PCA-S
srun gmx rdf -f traj_npt_main_nojump_100ps.trr -s py_tprfile -sf Pall.inp -o Pall.xvg
mv Pall.xvg ${rdfout_dir}
mv Pall.inp ${rdfout_dir}
wait

printf "End of RDF calculations.."
