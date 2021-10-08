# To run different initial conditions for polydisperse/monodisperse
# lignin for computing glass transition temperature
# Version: Apr-15-2021
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
from aux_tg_inpgen import * # function definitions
#------------------------------------------------------------------

# Version Info and # of input args for parsing thermostat coeff
print("Generating GROMACS run-time inputs")
print("Version: Apr-15-2021")
if len(sys.argv) == 1:
    coeff_fyle = 'None' 
elif len(sys.argv) == 2:
    coeff_fyle = str(sys.argv[1])
else:
    print('Unknown number of arguments: ', len(sys.argv),\
          str(sys.argv))
    exit()
#------------------------------------------------------------------

# Input Data
run_all   = 1 # 1-copy files and run, 0-NO run (copies files)
gpu       = 0 # 0-no gpu, 1 - gpu
inp_type  = 'melts' # melts, solvents, cosolvents
biomass   = 'WT' # name of the biomass type
disp_arr  = [1.0]
run_arr   = [2,3] # number of independent runs for a given biomass
high_temp = 600 # Run at high temperature for relaxation
temp_min  = 405 # Minimum temperature
temp_max  = 431 # Maximum temperature (< max; add +1 to desired)
temp_dt   = 10  # Temperature dt
npoly_res = 22  # number of polymer residues
box_dim   = 15  # Initial box size
solv_name = 'None' # add this later
wat_name  = 'None' # add this later
o_sol_typ = 'None' # add this later
#------------------------------------------------------------------

# Directory Paths
main_dir  = os.getcwd() # current dir
gmx_dir   = '../src_gmx' # gmx file super directory
top_dir   = gmx_dir + '/solv_files/topol' # topology dir
cfg_dir   = gmx_dir + '/solv_files/initguess' # configuration dir
itp_dir   = gmx_dir + '/solv_files/prm_itp' # prm file dir
mdp_dir   = gmx_dir + '/' + 'mdp_files' # mdp file dir
sh_dir    = gmx_dir + '/' + 'sh_files' # sh file dir
scr_dir   = '/lustre/or-scratch/cades-bsd/v0e' # scratch dir

if not os.path.isdir(scr_dir):
    print("FATAL ERROR: ", scr_dir, " not found")
    exit("Check scratch directory path")
scr_dir  = scr_dir + '/Glassy_lignin'
if not os.path.isdir(scr_dir):
    os.mkdir(scr_dir)
#------------------------------------------------------------------

# Required GMX/sh and default gro/top files
mdp_fyles  = ['minim_pyinp.mdp','nvt_pyinp.mdp','nvt_high_pyinp.mdp',\
              'npt_berendsen_pyinp.mdp','npt_main_pyinp.mdp']
if gpu:
    sh_md_fyle = 'run_md_gpu_pyinp.sh'
    sh_pp_fyle = 'run_preprocess_gpu_pyinp.sh'
else:
    sh_md_fyle = 'run_md_pyinp.sh'
    sh_pp_fyle = 'run_preprocess_pyinp.sh'
    
sh_rep_fyl = ['repeat_all.sh','repeat_md.sh']
def_inicon = 'initconf.gro'
#------------------------------------------------------------------

# Main code
for disp_val in range(len(disp_arr)): # loop in polydisperse array

    # Make directories
    head_dir = scr_dir + '/' + inp_type
    if not os.path.isdir(head_dir):
        print(head_dir, " does not exist")
        print("ERR: Create path")
        continue

    poly_dir = head_dir + '/' + biomass
    if not os.path.isdir(poly_dir):
        print(poly_dir, " does not exist")
        print("ERR: Create path")
        continue


    poly_dir = poly_dir + '/pdi_' + str(disp_arr[disp_val])
    if not os.path.isdir(poly_dir):
        print(poly_dir, " does not exist")
        print("ERR: Create path")
        continue

    for casenum in range(len(run_arr)): # loop in runarr

        rundir = poly_dir + '/run_' + str(run_arr[casenum])
        if not os.path.isdir(rundir):
            print(rundir, " does not exist")
            print('ERR: Create path and input top/pdb files first')
            continue

        
        # set main working directory and check input files are present
        workdir1 = set_working_dir(rundir,inp_type,o_sol_typ)

        # Loop over required temperature range
        for curr_temp in range(temp_min,temp_max,temp_dt):           

            poly_conffile,poly_topfile=check_inp_files(workdir1,'None')

            show_initial_log(biomass,o_sol_typ,disp_arr[disp_val],\
                             run_arr[casenum],curr_temp)

            temp_dir = workdir1 + '/T_'+str(curr_temp)
            if not os.path.isdir(temp_dir):
                os.mkdir(temp_dir)            

            # Set thermostat/top variables (change if there is a temp sweep)
            Tetau_nvt,Tetau_highnvt,Tetau_berend,Tetau_parrah,Prtau_berend,\
                Prtau_parrah,ref_pres,melt_topname=couple_coeff(inp_type,\
                                                                coeff_fyle)

            # Find tc_groups
            tc_grp,tc_typ=create_tcgrps(inp_type,npoly_res,solv_name,wat_name)

            # Copy and edit mdp files (all temperatures are same)
            check_cpy_mdp_files(mdp_dir,temp_dir,mdp_fyles,inp_type,Tetau_nvt\
                                ,Tetau_highnvt,Tetau_berend,Tetau_parrah,\
                                Prtau_berend,Prtau_parrah,curr_temp,high_temp,\
                                ref_pres,tc_grp,tc_typ,main_dir,coeff_fyle)
            
            # Check for tpr files. Re-edit run_md.sh irrespective of
            # whether tpr files are found. Edit run_pp.sh iff cont_run = 0
            cont_run = cpy_sh_files(sh_dir,temp_dir,sh_pp_fyle,sh_md_fyle)

            # Write index groups
            indx_fyle = create_indxgrps(temp_dir,inp_type,npoly_res,
                                        solv_name,wat_name)
            ff_dir = create_ff_dir(temp_dir) # empty dir for melts
 
            # Copy/edit top/conf/itp files for solvents/cosolvents
            if inp_type == 'solvents' or inp_type == 'cosolvents':
                if cont_run == 0:
                    # sol_top/sol_itp/sol_cfg are arrays
                    sol_top,sol_itp,sol_cfg=cpy_solv_files(top_dir,cfg_dir,\
                                                           itp_dir,ff_dir,\
                                                           inp_type,o_sol_typ,\
                                                           wat_type)
                    # check whether or not system is rerunning
                    poly_topedit = edit_main_top_file(poly_topfile,ff_dir,\
                                                      sol_top,sol_itp,temp_dir)
                    
                else: #already present
                    poly_topedit = poly_topfile

            else: # no need to change for melts
                poly_topedit = poly_topfile

            # Copy all shell script for running
            for fsh_name in sh_rep_fyl:
                #if not os.path.exists(temp_dir + '/' + fsh_name):
                gencpy(sh_dir,temp_dir,fsh_name) #copy always

            # Change to working directory
            os.chdir(temp_dir)

            # Make an outdir for writing job files from summit
            sum_out_dir = temp_dir + '/outdir'
            if not os.path.isdir(sum_out_dir):
                os.mkdir(sum_out_dir)

            # Edit pre-processing shell script files iff cont_run = 0
            if cont_run == 0: 
                # Copy polyconf/polytop for this casenum to sub-high_T dir    
                if curr_temp != high_temp: #for sub-high_T systems
                    new_conffile = find_hightemp_conf(workdir1,high_temp)
                    if new_conffile != 'none':
                        poly_conffile = new_conffile
                    
                #copy main polyconf file for high temperature and file from high temp sys for sub-high_T
                cpy_polyfiles_to_tempdir(poly_conffile,poly_topfile,workdir1,\
                                         temp_dir)
      
                if inp_type == 'melts':
                    edit_pp_files(biomass,inp_type,poly_conffile,0,0,\
                                  poly_topedit,'None','None',sh_pp_fyle\
                                  ,ff_dir,'None',0,indx_fyle,curr_temp)
                elif inp_type == 'solvents':
                    edit_pp_files(biomass,inp_type,poly_conffile,n_orgsolv\
                                  ,0,poly_topedit,o_sol_typ,'None',\
                                  sh_pp_fyle,ff_dir,sol_cfg,box_dim,\
                                  indx_fyle,curr_temp)
                elif inp_type == 'cosolvents':
                    edit_pp_files(biomass,inp_type,poly_conffile,n_orgsolv\
                                  ,nwater,poly_topedit,o_sol_typ,wat_type,\
                                  sh_pp_fyle,ff_dir,sol_cfg,box_dim,\
                                  indx_fyle,curr_temp)
            else: # check for initconf.gro
                if os.path.exists(def_inicon):
                    poly_conffile = def_inicon
                else:
                    print('ERR: tpr files found, but no', def_inicon, 'found')
                    continue

            # Edit run_md shell script files always
            edit_md_files(biomass,inp_type,poly_conffile,poly_topedit,\
                          o_sol_typ,sh_md_fyle,indx_fyle,curr_temp)

            # Submit preprocess job
            if run_all:
                if cont_run == 0:
                    print("Running all: repeat_all.sh")
                    subprocess.call(["chmod", "777", "repeat_all.sh"])
                    subprocess.call(["./repeat_all.sh"])
                else:
                    print("Running just MD part: repeat_md.sh")
                    subprocess.call(["chmod", "777", "repeat_md.sh"])
                    subprocess.call(["./repeat_md.sh"])

            print("Completed T = ", curr_temp)
            os.chdir(main_dir) #main dir

    # Write end of loop and return directory handle to main directory
    print( "End of run number: ", run_arr[casenum])
    os.chdir(main_dir) #main dir
