#!/bin/bash

#SBATCH -A bsd
#SBATCH -p burst
#SBATCH -t 0-23:00:00
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

#--------- Search for 100 ps interval file. If not present, create it ----------
if ! test -f "traj_npt_main_nojump_100ps.trr"; then
    if ! test -f "npt_main.tpr"; then
	printf "tpr/trr files not found"
	exit 1
    else
	printf "0" | srun gmx trjconv -s npt_main.tpr -f traj_npt_main.trr -dt 100 -pbc nojump -o traj_npt_main_nojump_100ps.trr
	wait
    fi
fi

#**********************Begin Rg/Eigenvalue calculation***************************
if py_rgflag ; then
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

fi
#----------------------End Rg/Eigenvalue calculation ---------------------------

#**********************Begin RDF calculation************************************
if py_rdfflag ; then
    # Inputs
    rdfout="rdf_nptmain"; rdfout_dir="all_rdfs"
    nchains=py_nchains
    printf "Computing RDF of chains"
    mkdir -p ${rdfout_dir}
       
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

fi
#----------------------End RDF calculation -------------------------------------

#**********************Begin HB calculation*************************************
if py_hbflag ; then
    # Inputs
    hbout="hb_nptmain"; hbout_dir="all_hbs"
    nchains=py_nchains
    printf "Computing hydrogen bonding"
    mkdir -p ${hbout_dir}

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

fi
#----------------------End HB --------------------------------------------------

#**********************Begin Segmental Rg calculation***************************
if py_segrgflag ; then
    # Inputs
    rgout="traj_npt_main.trr"; allresultdir="all_results"
    nchains=py_nchains

    # Search for compRg.psf in the superdirectory
    if ! test -f "../init_files/compRg.psf"; then
	printf "psf file not found"
	exit 1
    fi
    
    
    # Compute segmental Rg of chains
    printf "Computing segmental Rg of chains"
    
    if [ ! -d ${allresultdir} ]; then
	mkdir -p ${allresultdir}
    fi
    wait
    
    # Segmental Rg computation command using tcl script
    vmd -dispdev text -e calc_seg_rg.tcl
    
    mv rg_allsegs.dat ${allresultdir}
    mv RgvsN.dat ${allresultdir}

    printf "End of segmental Rg calculations.."
fi
#----------------------End Segmental Rg calculation ----------------------------

#**********************Begin MSD calculation************************************
if py_msdflag ; then

    # Inputs
    msdout="msd_nptmain"; allresultdir="all_results"
    nchains=py_nchains

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
fi
