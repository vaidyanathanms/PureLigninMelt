# To generate initial directories and fill with initial files for
# lignin systems. Will run LigninBuilder and generate initial
# pdb/top/psf files

# Version: Mar-14-2021
#------------------------------------------------------------------

# Import modules
import os
import sys
import numpy
import re
import shutil
import glob
import math
import fileinput
from supp_initdirs import *
#------------------------------------------------------------------

# Version Info and # of input args for parsing Lignin Builder
print("Generating GROMACS run-time inputs")
print("Version: Mar-03-2021")
if len(sys.argv) == 2:
    inp_fyle = str(sys.argv[1])
else:
    print('Unknown number of arguments: ', len(sys.argv),\
          str(sys.argv))
    exit()
#------------------------------------------------------------------

# Input Data
inp_type  = 'melts' # melts, solvents, cosolvents
biomass   = 'WT' # name of the biomass type
disp_arr  = [1.8]
run_arr   = [2] # number of independent runs for a given biomass
npoly_res = 22  # number of polymer residues
num_chain = 20  # number of polymer chains
solv_name = 'None' # add this later
wat_name  = 'None' # add this later
o_sol_typ = 'None' # add this later
#------------------------------------------------------------------

# Directory Paths
main_dir  = os.getcwd() # current dir
scr_dir   = '/lustre/or-scratch/cades-bsd/v0e' # scratch dir

if not os.path.isdir(scr_dir):
    print("FATAL ERROR: ", scr_dir, " not found")
    exit("Check scratch directory path")

scr_dir  = scr_dir + '/Glassy_lignin'
if not os.path.isdir(scr_dir):
    os.mkdir(scr_dir)
#------------------------------------------------------------------

# Required GMX/sh Files
lig_fyles = ['make_genpsf.py','genconf.py','findmissingterms.py']
#------------------------------------------------------------------

for disp_cnt in range(len(disp_arr)): # loop in polydisperse array

    print( "Initializing for PDI: " + str(disp_arr[disp_cnt]))
    # Make directories
    head_dir = scr_dir + '/' + inp_type
    if not os.path.isdir(head_dir):
        os.mkdir(head_dir)

    poly_dir = head_dir + '/' + biomass
    if not os.path.isdir(poly_dir):
        os.mkdir(poly_dir)


    poly_dir = poly_dir + '/pdi_' + str(disp_arr[disp_cnt])
    if not os.path.isdir(poly_dir):
        os.mkdir(poly_dir)

    for run_id in range(len(run_arr)): # loop in runarr

        print( "Initializing for run_" + str(run_arr[run_id]))
        outdir = poly_dir + '/run_' + str(run_arr[run_id])
        if not os.path.isdir(outdir):
            os.mkdir(outdir)


        # Retrieve case number for moving into the new directory
        casenum,clrdir = retrieve_case_num(inp_fyle)
        if clrdir == 1: #clear directory if they exist
            init_dir = main_dir + '/casenum_' + str(casenum)
            if os.path.isdir(init_dir):
                shutil.rmtree(init_dir)

        # Edit generic python input file
        edited_inpfyle = editinpfyle(main_dir,inp_fyle,\
                                     disp_arr[disp_cnt],run_arr[run_id],\
                                     npoly_res,num_chain)
        # Run genconf.py
        run_genconf(edited_inpfyle,lig_fyles,casenum)
        # Check casenum directory is made
        init_dir = main_dir + '/casenum_' + str(casenum)
        if not os.path.isdir(init_dir):
            raise RuntimeError(init_dir + " does not exist!")

        os.chdir(init_dir)
        # Run all tcl steps from init_dir
        run_all_steps(init_dir)

        if inp_type == 'melts':
            findir = set_working_dir(outdir,inp_type,o_sol_typ)
            supp_findir = findir + '/init_files'
            if not os.path.isdir(supp_findir):
                os.mkdir(supp_findir)

            # Copy all pdb/topology files
            cpy_pdb_top(init_dir,findir,biomass)
            cpy_supp_files(init_dir,supp_findir,disp_arr[disp_cnt])

        else:
            # Copy same biomass to all solvent directories
            for solinc in range(len(o_sol_arr)):
                o_sol_typ = o_sol_arr[solinc]
                findir = set_working_dir(outdir,inp_type,o_sol_typ)
                supp_findir = findir + '/init_files'
                if not os.path.isdir(supp_findir):
                    os.mkdir(supp_findir)

                # Copy all pdb/topology files
                cpy_pdb_top(init_dir,findir,biomass)
                cpy_supp_files(init_dir,supp_findir)

        # End details and return directory handle
        print("End making files for run_" + str(run_arr[run_id]))
        os.chdir(main_dir)
