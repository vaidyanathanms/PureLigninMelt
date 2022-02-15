# Generate GROMACS inputs for analysis
# Ver: May-11-2021
#------------------------------------------------------------------

# Import modules
import os
import sys
import numpy
import re
import shutil
import glob
import math
import subprocess
from aux_gmx_functions import * # function definitions
#------------------------------------------------------------------

# Version Info and parse input conditions
print("Generating GROMACS analysis inputs")
print("Version: May-11-2021")
#------------------------------------------------------------------

# Input Keys
rg_calc    = 0 # Calculate rg
seg_rgcalc = 0 # Calculate segmental rg
msd_calc   = 0 # Calculate msd
rdf_calc   = 0 # Calculate rdf
shape_calc = 0 # Calculate shape factor
tg_calc    = 1 # Calculate densities for Tg
#------------------------------------------------------------------

# Input Data
run_all   = 1 # 1-copy files and run, 0-NO run (copies files)
inp_type  = 'melts' # melts, solvents, cosolvents
biomass   = 'WT' # name of the biomass type
disp_arr  = [1.8]# dispersity values
run_arr   = [2,3]  # run number for a given dispersity
temp_min  = 300  # Minimum temperature
temp_max  = 321  # Maximum temperature (< max; add +1 to desired)
temp_dt   = 20   # Temperature dt
nchains   = 20   # Number of chains - cross check from conf file
solv_name = 'None' # add this later
wat_name  = 'None' # add this later
o_sol_typ = 'None' # add this later
#------------------------------------------------------------------

# Directory Paths
main_dir  = os.getcwd() # current dir
gmx_dir   = '/home/v0e/allcodes/files_lignin/pure_melts/myruns/src_gmx' # gmx dir
tcl_dir   = '/home/v0e/allcodes/files_lignin/pure_melts/myruns/src_tcl' # tcl dir
sh_dir    = gmx_dir + '/' + 'sh_files' # sh file dir
scr_dir   = '/lustre/or-scratch/cades-bsd/v0e' # scratch dir
#------------------------------------------------------------------

if not os.path.isdir(scr_dir):
    print("FATAL ERROR: ", scr_dir, " not found")
    exit("Check scratch directory path")
scr_dir  = scr_dir + '/Glassy_lignin'
if not os.path.isdir(scr_dir):
    os.mkdir(scr_dir)
#------------------------------------------------------------------

# Main code
for disp_val in range(len(disp_arr)): # loop in polydisperse array

    # Make directories
    head_dir = scr_dir + '/' + inp_type
    if not os.path.isdir(head_dir):
        print(head_dir, " does not exist")
        continue

    poly_dir = head_dir + '/' + biomass
    if not os.path.isdir(poly_dir):
        print(poly_dir, " does not exist")
        continue

    poly_dir = poly_dir + '/pdi_' + str(disp_arr[disp_val])
    if not os.path.isdir(poly_dir):
        print(poly_dir, " does not exist")
        continue

    for casenum in range(len(run_arr)): # loop in runarr

        rundir = poly_dir + '/run_' + str(run_arr[casenum])
        if not os.path.isdir(rundir):
            print(rundir, " does not exist")
            continue
        
        # set main working directory and check input files are present
        workdir1 = set_working_dir(rundir,inp_type,o_sol_typ)

        # Tg calculation is for entire temperature range
        if tg_calc:
            if not os.path.isdir(workdir1 + '/outdir'):
                os.mkdir(outdir)
            tlist = [temp_min,temp_max,temp_dt]
            set_jobfile(nchains,sh_dir,workdir1,'comp_dens_pyinp.sh',\
                        'None','None','None',tlist,'rho')

        # Loop over required temperature range
        for curr_temp in range(temp_min,temp_max,temp_dt): 

            temp_dir = workdir1 + '/T_'+str(curr_temp)
            if not os.path.isdir(temp_dir):
                print(temp_dir, " does not exist")
                continue

            print('Analyzing: ', biomass,inp_type,'run_'+\
                  str(run_arr[casenum]), 'T_'+str(curr_temp))

            conffile,topfile=check_inp_files(temp_dir,'None')
            mon_list,at_list = count_nchains(temp_dir,conffile,\
                                             nchains)
            tprfile, trajfile = find_tpr_trr_files(temp_dir)
            if tprfile == -1 or trajfile == -1:
                print("ERROR: Did not find tpr/trr file...")
                continue

            create_anagrps_inp(temp_dir,mon_list,at_list,nchains)

            # Copy and set scripts for running
            if rg_calc:
                set_jobfile(nchains,sh_dir,temp_dir,'rgcomp_pyinp.sh',\
                            trajfile,tprfile,conffile,curr_temp,'rg')
            if seg_rgcalc:
                set_jobfile(nchains,sh_dir,temp_dir,'rgcomp_pyinp.sh',\
                            trajfile,tprfile,conffile,curr_temp,'segrg')
            if msd_calc:
                set_jobfile(nchains,sh_dir,temp_dir,'msdcomp_pyinp.sh',\
                            trajfile,tprfile,conffile,curr_temp,'msd')
            if rdf_calc:
                set_jobfile(nchains,sh_dir,temp_dir,'rdfcomp_pyinp.sh',\
                            trajfile,tprfile,conffile,curr_temp,'rdf')
            if shape_calc:
                set_jobfile(nchains,sh_dir,temp_dir,'shapecomp_pyinp.sh',\
                            trajfile,tprfile,conffile,curr_temp,'shape')

            os.chdir(temp_dir) # Change to temperature directory
            if run_all:
                if tg_calc:
                    os.chdir(workdir1)
                    print("Running Tg calculation")
                    subprocess.call(["sbatch","comp_dens.sh"])
                    os.chdir(temp_dir)
                if rg_calc:
                    print("Running Rg calculation")
                    subprocess.call(["sbatch","rgcomp.sh"])
                if seg_rgcalc:
                    print("Running segmental Rg calculation")
                    gencpy(tcl_dir,temp_dir,'calc_seg_rg.tcl')
                    subprocess.call(["sbatch","segment_rgcomp.sh"])
                if msd_calc:
                    print("Running MSD calculation")
                    subprocess.call(["sbatch","msdcomp.sh"])
                if rdf_calc:
                    print("Running RDF calculation")
                    subprocess.call(["sbatch","rdfcomp.sh"])
                if shape_calc:
                    print("Running shape factor calculation")
                    subprocess.call(["sbatch","shapecomp.sh"])
            print("Completed T = ", curr_temp)
            os.chdir(main_dir) #main dir
