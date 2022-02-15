# Version: May-22-2020
# To consolidate outputs and restart files

import numpy
import os
import shutil
import subprocess
import sys
import glob

#---------------Import Functions-----------------------------------

from subprocess import call
def my_cpy_generic(srcdir,destdir,inpfylname,destfylname):
    src_fyl  = srcdir  + '/' + inpfylname
    dest_fyl = destdir + '/' + destfylname
    shutil.copy2(src_fyl, dest_fyl)

#---------input flags------------------------------------------
analysis_type  = 2 #1-polydisp_pe,2-newparams,3-mwchange,#4-oldsize
rest_job_flag  = 0 #copy restart/job files
fyl_flag       = 1 #copy output files

#---------input details----------------------------------------
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

#--------file_lists--------------------------------------------
#Give prefix for files to be copied followed by *
fyl_list     = ['adsfrac*','tether*','chainadsval*','log*',\
                'dens*','chdens*','PErdf*','chgrpdens*','grpdens*'\
                ,'polydens*','anainp*']
restart_list = ['restart*','archival*','job*']

#---------directory info---------------------------------------
main_dir  = os.getcwd() # current dir
gmx_dir   = '/home/v0e/allcodes/files_lignin/pure_melts/myruns/src_gmx' # gmx dir
tcl_dir   = '/home/v0e/allcodes/files_lignin/pure_melts/myruns/src_tcl' # tcl dir
sh_dir    = gmx_dir + '/' + 'sh_files' # sh file dir
scr_dir   = '/lustre/or-scratch/cades-bsd/v0e' # scratch dir

if not os.path.isdir(scr_dir):
    print(scr_dir, "does not exist")
    sys.exit()
#--------------Main Analysis------------------------------------

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
        
        #------------------Make global analysis files output directory-----
        if fyl_flag == 1:

            out_dir = 'results_PDI_' + str(disp_arr[disp_val])         
            anafyl_maindir = poly_dir + '/' + out_dir
            if not os.path.isdir(anafyl_maindir):
                print("making", anafyl_maindir)
                os.mkdir(anafyl_maindir)
                
            anafyl_casedir = anafyl_maindir + '/run_' + str(run_arr[casenum])
            if not os.path.isdir(anafyl_casedir):
                print("making", anafyl_casedir)
                os.mkdir(anafyl_casedir)

        if rest_job_flag == 1:

            rest_dirname = 'rest_PDI_' + str(disp_arr[disp_val])
            rest_maindir = poly_dir + '/' + rest_dirname
            if not os.path.isdir(rest_maindir):
                os.mkdir(restart_maindir)

            rest_casedir = rest_maindir + '/run_' + str(run_arr[casenum])
            if not os.path.isdir(rest_casedir):
                print("making", restart_casedir)
                os.mkdir(restart_casedir)

        for curr_temp in range(temp_min,temp_max,temp_dt): 

            temp_dir = workdir1 + '/T_'+str(curr_temp)
            if not os.path.isdir(temp_dir):
                print(temp_dir, " does not exist")
                continue

            os.chdir(temp_dir)
            destdir = os.getcwd()
            print('Copying: ', biomass,inp_type,'run_'+\
                  str(run_arr[casenum]), 'T_'+str(curr_temp))


            #----Make temp directory in pdi_directory----------
            if fyl_flag == 1:
                anafyl_dir = anafyl_casedir + '/T_' + str(curr_temp))
                if not os.path.isdir(anafyl_dir):
                    os.mkdir(anafyl_dir)
            if rest_job_flag == 1:
                restart_dir = rest_casedir + '/T_' + str(curr_temp))
                if not os.path.isdir(restart_dir):
                    os.mkdir(restart_dir)

            #------All copying/manipulations--------------------------
            result_dir = 'results_'+ str(free_chains[ifree]) + '_' \
                         + str(pdi_free) + '_' + str(ncases_pdi[casenum])
            workdir_results = workdir_subpdi + '/' + result_dir
            if not os.path.isdir(workdir_results):
                print(workdir_results, " does not exist")
                continue

            print("Copying output files from", workdir_results)

            if fyl_flag == 1:

                for fylcnt in range(len(fyl_list)):

                    #search in results_directory
                    os.chdir(workdir_results)
                    destdir = os.getcwd()
                    list_of_files = glob.glob(fyl_list[fylcnt])

                    if list_of_files == []: #if list is empty in results dir

                        #search in previous directory
                        os.chdir(workdir_subpdi)
                        destdir = os.getcwd()
                        list_of_files = glob.glob(fyl_list[fylcnt])

                        if list_of_files == []: #if list is empty here too
                            print("Did not find files of type", \
                                  fyl_list[fylcnt])
                            continue

                    print( "Copying files of type: ", fyl_list[fylcnt])

                    for filenum in range(len(list_of_files)):

                        fylname = list_of_files[filenum]
                        anafylname = fylname
                        my_cpy_generic(destdir,anafyl_dir,fylname,anafylname)

            if rest_job_flag == 1:

                for fylcnt in range(len(restart_list)):

                    #search only in main directory
                    os.chdir(workdir_subpdi)
                    destdir = os.getcwd()
                    list_of_files = glob.glob(restart_list[fylcnt])

                    if list_of_files == []: #
                        print("Did not find files of type", \
                              restart_list[fylcnt])
                        continue

                    print( "Copying files of type: ", restart_list[fylcnt])
                    for filenum in range(len(list_of_files)):

                        fylname = list_of_files[filenum]
                        anafylname = fylname
                        my_cpy_generic(destdir,restart_dir,fylname,anafylname)# Aux
