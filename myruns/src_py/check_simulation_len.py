# Aux python codes for generating GROMACS inputs for analysis
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
import time
from general_functions import * # function definitions
#------------------------------------------------------------------

# Version Info and parse input conditions
print("Generating GROMACS analysis inputs")
print("Version: May-11-2021")
#------------------------------------------------------------------

# Input Data
expt      = 0 # Experimental distribution - 1
inp_type  = 'melts' # melts, solvents, cosolvents
biomass   = 'WT' # name of the biomass type
disp_arr  = [1.0,1.8,3.0,3.7]# dispersity values
run_arr   = [1,2,3,4]  # run number for a given dispersity
temp_min  = 300  # Minimum temperature
temp_max  = 501  # Maximum temperature (< max; add +1 to desired)
temp_dt   = 10   # Temperature dt
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
fout = open("run_stats.dat", "a") #append to stats file
datetime = time.strftime("%m/%d/%Y, %H:%M:%S")
fout.write("**Date and Time** %s\n" %(datetime))
fout.write("%s\t%s\t%s\t%s\t%s\n" \
           %("PDI","Run","Temp","Filename","Time"))

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
    
    if expt:
        poly_dir = poly_dir + '/expts'
    else:
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
                fout.write('%g\t%d\t%g\t%s\t%d\n' \
                           %(disp_arr[disp_val],run_arr[casenum],\
                             curr_temp,'None',-1))
                continue
                
            outfile_dir = temp_dir + '/trajfiles'
            if not os.path.isdir(outfile_dir):
                print(temp_dir, " does not exist")
                fout.write('%g\t%d\t%g\t%s\t%d\n' \
                           %(disp_arr[disp_val],run_arr[casenum],\
                             curr_temp,'None',-1))
                continue

            logout = find_latest_file(outfile_dir,'.log')


            if logout == 'None': # if no log file is found
                fout.write('%g\t%d\t%g\t%s\t%d\n' \
                           %(disp_arr[disp_val],run_arr[casenum],\
                             curr_temp,'None',-1))
                continue
            
            for line in read_reverse_order(logout):
                if "Statistics" in line:
                    if expt:
                        fout.write('%s\t%d\t%g\t%s\t%s\n' \
                                   %('Expt' + str(disp_arr[disp_val])\
                                     ,run_arr[casenum],curr_temp,\
                                     extract_file_name(logout),\
                                     line.strip().split()[2]))
                        break
                    else:
                         fout.write('%g\t%d\t%g\t%s\t%s\n' \
                                    %(disp_arr[disp_val],\
                                      run_arr[casenum],curr_temp,\
                                      extract_file_name(logout),\
                                      line.strip().split()[2]))
                         break

        
fout.close()
            


