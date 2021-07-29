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
if len(sys.argv) == 1:
    coeff_fyle = 'None' 
elif len(sys.argv) == 2:
    coeff_fyle = str(sys.argv[1])
else:
    print('Unknown number of arguments: ', len(sys.argv),\
          str(sys.argv))
    exit()
#------------------------------------------------------------------

# Input Keys
rg_calc   = 1 # Calculate rg
msd_calc  = 0 # Calculate msd
rdf_calc  = 0 # Calculate rdf

# Input Data
run_all   = 1 # 1-copy files and run, 0-NO run (copies files)
inp_type  = 'melts' # melts, solvents, cosolvents
biomass   = 'WT' # name of the biomass type
disp_arr  = [3.0]
run_arr   = [6] # number of independent runs for a given biomass
temp_min  = 360 # Minimum temperature
temp_max  = 501 # Maximum temperature (< max; add +1 to desired)
temp_dt   = 20  # Temperature dt
npoly_res = 22  # number of polymer residues
nchains   = 20  # Number of chains - cross check from conf file
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
#------------------------------------------------------------------

if not os.path.isdir(scr_dir):
    print("FATAL ERROR: ", scr_dir, " not found")
    exit("Check scratch directory path")
scr_dir  = scr_dir + '/Glassy_lignin'
if not os.path.isdir(scr_dir):
    os.mkdir(scr_dir)
#------------------------------------------------------------------

#v Required GMX/sh and default gro/top files
def_inicon = 'initconf.gro'
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
                continue

            create_anagrps_inp(temp_dir,mon_list,at_list,nchains)
            run_analysis(nchains,rg_calc,msd_calc,rdf_calc,sh_dir,\
                         temp_dir,trajfile,tprfile,conffile,curr_temp)

            os.chdir(temp_dir)
            if run_all:
                if rg_calc:
                    print("Running Rg calculation")
                    subprocess.call(["sbatch","rgcomp.sh"])
                if msd_calc:
                    print("Running MSD calculation")
                    subprocess.call(["sbatch","msdcomp.sh"])
                if rdf_calc:
                    print("Running RDF calculation")
                    subprocess.call(["sbatch","rdfcomp.sh"])

            print("Completed T = ", curr_temp)
            os.chdir(main_dir) #main dir
