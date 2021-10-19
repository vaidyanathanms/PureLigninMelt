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
from aux_gmx_functions import gencpy 
from aux_gmx_functions import set_working_dir
#------------------------------------------------------------------

# Version Info and parse input conditions
print("Copying outputs")
print("Version: Oct-04-2021")
#------------------------------------------------------------------

# Input Keys
rg_calc    = 1 # Calculate rg
seg_rgcalc = 1 # Calculate segmental rg
msd_calc   = 1 # Calculate msd
rdf_calc   = 0 # Calculate rdf
#------------------------------------------------------------------

# Input Data
inp_type  = 'melts' # melts, solvents, cosolvents
biomass   = 'WT' # name of the biomass type
disp_arr  = [1.8]# dispersity values
run_arr   = [1]  # run number for a given dispersity
temp_min  = 300  # Minimum temperature
temp_max  = 501  # Maximum temperature (< max; add +1 to desired)
temp_dt   = 20   # Temperature dt
solv_name = 'None' # add this later
wat_name  = 'None' # add this later
o_sol_typ = 'None' # add this later
#------------------------------------------------------------------

# Directory Paths
main_dir  = os.getcwd() # current dir
gmx_dir   = '/home/v0e/allcodes/files_lignin/pure_melts/myruns/src_gmx' # gmx dir
tcl_dir   = '/home/v0e/allcodes/files_lignin/pure_melts/myruns/src_tcl' # tcl dir
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

    headout2_dir  = poly_dir + '/all_outputs'
    if not os.path.isdir(headout2_dir):
        os.mkdir(headout2_dir)

    poly_dir = poly_dir + '/pdi_' + str(disp_arr[disp_val])
    if not os.path.isdir(poly_dir):
        print(poly_dir, " does not exist")
        continue

    headout_dir = headout2_dir + '/pdi_' + str(disp_arr[disp_val])
    if not os.path.isdir(headout_dir):
        os.mkdir(headout_dir)

    for casenum in range(len(run_arr)): # loop in runarr

        rundir = poly_dir + '/run_' + str(run_arr[casenum])
        if not os.path.isdir(rundir):
            print(rundir, " does not exist")
            continue

        runout_dir = headout_dir + '/run_' + str(run_arr[casenum])
        if not os.path.isdir(runout_dir):
            os.mkdir(runout_dir)
        
        #set main working directory and check input files are present
        workdir1 = set_working_dir(rundir,inp_type,o_sol_typ)
        if not os.path.isdir(workdir1):
            print(workdir1, " does not exist")
            continue

        # Loop over required temperature range
        for curr_temp in range(temp_min,temp_max,temp_dt): 

            temp_dir = workdir1 + '/T_'+str(curr_temp)
            if not os.path.isdir(temp_dir):
                print(temp_dir, " does not exist")
                continue

            print('Copying data from: ', biomass,inp_type,'run_'+\
                  str(run_arr[casenum]), 'T_'+str(curr_temp))
            
            if rg_calc or seg_rgcalc:
                print("Copying Rg files")

                rg_headout  = runout_dir + '/rganalysis'
                if not os.path.isdir(rg_headout):
                    os.mkdir(rg_headout)

                rg_indir = temp_dir + '/rganalysis'
                if not os.path.isdir(rg_indir):
                    print(rg_indir, "does not exist")
                    continue

                rg_outdir = rg_headout + '/T_' + str(curr_temp)
                if not os.path.isdir(rg_outdir):
                    os.mkdir(rg_outdir)

                list_rg = glob.glob(rg_indir + '/*')
                if list_rg == []:
                    print('No files found inside ', rg_indir)
                    continue
                
                for f in list_rg:
                    fname = f.split('/')[len(f.split('/'))-1]
                    gencpy(rg_indir,rg_outdir,fname)

            if msd_calc:

                print("Copying MSD files")

                msd_headout  = runout_dir + '/msdanalysis'
                if not os.path.isdir(msd_headout):
                    os.mkdir(msd_headout)

                msd_indir = temp_dir + '/msdanalysis'
                if not os.path.isdir(msd_indir):
                    print(msd_indir, "does not exist")
                    continue

                msd_outdir = msd_headout + '/T_' + str(curr_temp)
                if not os.path.isdir(msd_outdir):
                    os.mkdir(msd_outdir)

                list_msd = glob.glob(msd_indir + '/*')
                if list_msd == []:
                    print('No files found inside ', msd_indir)
                    continue
                
                for f in list_msd:
                    fname = f.split('/')[len(f.split('/'))-1]
                    gencpy(msd_indir,msd_outdir,fname)

                
            print("Completed T = ", curr_temp)
            os.chdir(main_dir) #main dir
