# Version: May-22-2020
# To consolidate outputs and restart files

#---------------Import Functions-----------------------------------
import numpy
import os
import shutil
import subprocess
import sys
import glob
from subprocess import call
from general_functions import *
#---------input flags------------------------------------------
rest_job_flag  = 1 #copy restart/job files
fyl_flag       = 1 #copy output files

#---------input details----------------------------------------
inp_type  = 'melts' # melts, solvents, cosolvents
biom_arr  = ['WT']  # name of the biomass type
pdi_arr   = [1.0,1.8]#,3.7,'expts'] # dispersity values
run_arr   = [1]#,2,3,4]         # run number for a given dispersity
temp_min  = 250     # Minimum temperature
temp_max  = 271     # Maximum temperature (< max; add +1 to desired)
temp_dt   = 10      # Temperature dt

#--------file_lists--------------------------------------------
#Give prefix for files to be copied with *
dir_list     = ['all_hbs','all_radgyr','all_results','all_rdfs']
fyl_list     = ['*.xvg','*.inp','*.dat','*.ndx']
restart_list = ['*.cpt','*.log','*.gro']

#---------directory info---------------------------------------
main_dir  = os.getcwd() # current dir
gmx_dir   = '/home/v0e/allcodes/files_lignin/pure_melts/myruns/src_gmx' # gmx dir
tcl_dir   = '/home/v0e/allcodes/files_lignin/pure_melts/myruns/src_tcl' # tcl dir
sh_dir    = gmx_dir + '/' + 'sh_files' # sh file dir
scr_dir   = '/lustre/or-scratch/cades-bsd/v0e/Glassy_lignin' # scratch dir

if not os.path.isdir(scr_dir):
    print(scr_dir, "does not exist")
    sys.exit()
#--------------Main Analysis------------------------------------

for bio_indx in range(len(biom_arr)): # loop in biomass
    
    biomass = biom_arr[bio_indx]
    headdir = scr_dir + '/' + inp_type + '/'+ biomass
    if not os.path.isdir(headdir):
        print("ERR: " + headdir + " not found")
        continue

    for pdi_val in pdi_arr:

        #--------------Create top-level output directories---------------
        allop_main = headdir + '/results_pdi_' + str(pdi_val)
        if not os.path.isdir(allop_main):
            os.mkdir(allop_main)

        supp_dir = scr_dir + '/' + inp_type + '/'+ biomass +\
                   '/restartfiles_pdi_' + str(pdi_val)
        if not os.path.isdir(supp_dir):
            os.mkdir(supp_dir)

      
        for casenum in range(len(run_arr)): # loop in runarr
            
            if pdi_val == 'expts':
                rundir = headdir + '/expts' + '/run_' + str(run_arr[casenum])
            else:
                rundir = headdir + '/pdi_' + str(pdi_val)+ '/run_'\
                         + str(run_arr[casenum])
            if not os.path.isdir(rundir):
                print(rundir + " does not exist!")
                continue
            #--------------Make global analysis files output directory-----
            if fyl_flag == 1:
                
                ana_casedir = allop_main + '/run_' + str(run_arr[casenum])
                if not os.path.isdir(ana_casedir):
                    print("making", ana_casedir)
                    os.mkdir(ana_casedir)
                    
            if rest_job_flag == 1:

                rest_casedir = supp_dir + '/run_' + str(run_arr[casenum])
                if not os.path.isdir(rest_casedir):
                    print("making", rest_casedir)
                    os.mkdir(rest_casedir)

            for tval in range(temp_min,temp_max,temp_dt): # loop in temp
                
                tdir = rundir + '/T_' + str(tval)

                print('Copying: ',biomass,inp_type,'run_'+str(run_arr[casenum]),\
                      'T_'+str(tval))

                #----Make temp directory in result directory---------------------
                if fyl_flag == 1: #output results
                    anafyl_dir = ana_casedir + '/T_' + str(tval)
                    if not os.path.isdir(anafyl_dir):
                        os.mkdir(anafyl_dir)
                if rest_job_flag == 1: #restart backup
                    restart_dir = rest_casedir + '/T_' + str(tval)
                    if not os.path.isdir(restart_dir):
                        os.mkdir(restart_dir)

                #------Copying output/restart files ------------------------------
                if fyl_flag:

                    for dirname in dir_list: #sub-directories in output
                        workdir_results = tdir + '/' + dirname
                        if not os.path.isdir(workdir_results):
                            print('Did not find the directory + ', workdir_results)
                            continue

                        qty_outdir = anafyl_dir + '/' + dirname
                        if not os.path.isdir(qty_outdir):
                            os.mkdir(qty_outdir)

                        print("Copying output files from", workdir_results)
                        for fylcnt in range(len(fyl_list)): #search for filetype
                            # search in subdirectory type
                            os.chdir(workdir_results)
                            desdir = os.getcwd()
                            list_of_files = glob.glob(fyl_list[fylcnt])

                            if list_of_files == []: #if list is empty

                                print("No filetype: ",fyl_list[fylcnt])
                                continue

                            print( "Copying files of type: ", fyl_list[fylcnt])

                            for filenum in range(len(list_of_files)):

                                fylname = list_of_files[filenum]
                                gencpy(desdir,qty_outdir,fylname)

                if rest_job_flag:

                    print("Copying restart files from", tdir)
                    for fylcnt in range(len(restart_list)):

                        #search only in main directory
                        os.chdir(tdir)
                        desdir = os.getcwd()
                        list_of_files = glob.glob(restart_list[fylcnt])

                        if list_of_files == []: #
                            print("Did not find files of type", \
                                  restart_list[fylcnt])
                            continue

                        print( "Copying files of type: ", restart_list[fylcnt])
                        for filenum in range(len(list_of_files)):

                            fylname = list_of_files[filenum]
                            gencpy(desdir,restart_dir,fylname)
