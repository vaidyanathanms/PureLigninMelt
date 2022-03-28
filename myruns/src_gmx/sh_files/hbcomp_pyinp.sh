#!/bin/bash

#SBATCH -A bsd
#SBATCH -p burst
#SBATCH -t 0-07:00:00
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

# Compute inter/intra HB
for (( chcnt_i = 0; chcnt_i < nchains-1; chcnt_i++ ))
do
    # First make corresponding index file
    srun gmx select -s py_tprfile -on hb${chcnt_i}.ndx -sf hb${chcnt_i}.inp
    wait
    # Now compute intra and inter for each chain
    echo "0 1" | srun gmx hbond -f traj_npt_main_nojump_100ps.trr -s py_tprfile -n hb${chcnt_i}.ndx -num ${hbout}_intra_${chcnt_i}.xvg
wait
    echo "0 2" | srun gmx hbond -f traj_npt_main_nojump_100ps.trr -s py_tprfile -n hb${chcnt_i}.ndx -num ${hbout}_inter_${chcnt_i}.xvg
wait
done

mv ${hbout}_*.xvg ${hbout_dir}
mv hb*.ndx ${hbout_dir}
mv hb*.inp ${hbout_dir}
wait

# Compute FA-FA/FA-PCA/FA-H/FA-G/FA-S
# Make index files
if ! test -f resinp.ndx; then
    if ! test -f "resinp.inp"; then
	echo "resinp.inp for residue indices not found!"
	exit 1
    else
	srun gmx select -s py_tprfile -on resinp.ndx -sf resinp.inp
	wait
    fi
fi

# Compute HBs
echo "0 0" | srun gmx hbond -f traj_npt_main_nojump_100ps.trr -s py_tprfile -n resinp.ndx -num hb_ff.xvg
wait
echo "0 1" | srun gmx hbond -f traj_npt_main_nojump_100ps.trr -s py_tprfile -n resinp.ndx -num hb_fp.xvg
wait
echo "0 2" | srun gmx hbond -f traj_npt_main_nojump_100ps.trr -s py_tprfile -n resinp.ndx -num hb_fh.xvg
wait
echo "0 3" | srun gmx hbond -f traj_npt_main_nojump_100ps.trr -s py_tprfile -n resinp.ndx -num hb_fg.xvg
wait
echo "0 4" | srun gmx hbond -f traj_npt_main_nojump_100ps.trr -s py_tprfile -n resinp.ndx -num hb_fs.xvg
wait
mv hb_f*.xvg ${hbout_dir}
mv resinp.inp ${hbout_dir}
mv resinp.ndx ${hbout_dir}
wait

printf "End of HB calculations.."
