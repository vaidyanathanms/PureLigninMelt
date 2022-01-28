# Supplementary python files for initialize_dirs_for_runs.py
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
import subprocess
#------------------------------------------------------------------

# General copy script
def gencpy2(dum_maindir,dum_destdir,fylname,fylname2):

    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        print('ERROR: ', srcfyl, 'not found')
        return

    desfyl = dum_destdir + '/' + fylname2
    shutil.copy2(srcfyl,desfyl)
#------------------------------------------------------------------

# General copy script
def gencpy(dum_maindir,dum_destdir,fylname):

    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        print('ERROR: ', srcfyl, 'not found')
        return

    desfyl = dum_destdir + '/' + fylname
    shutil.copy2(srcfyl,desfyl)
#------------------------------------------------------------------

# Set working directory
def set_working_dir(rundir,inp_type,solv_type = 'None'):
    if inp_type == 'solvents':
        workdir1 = rundir + '/' + solv_type
    elif inp_type == 'cosolvents':
        h_workdir = rundir + '/' + solv_type
        if not os.path.isdir(h_workdir):
            os.mkdir(h_workdir)
        workdir1 = h_workdir + '/' + solv_type + '_water'
    else:
        workdir1 = rundir #even for inp_type = melts

    if not os.path.isdir(workdir1):
        os.mkdir(workdir1)
    
    return workdir1
#------------------------------------------------------------------

# Retrieve case number
def retrieve_case_num(coeff_fyle):
    casenum = 1 # default
    clear_all = 1 # clear_all files
    with open(coeff_fyle) as farg:
        for line in farg:
            line = line.rstrip('\n')
            if line.startswith('#'):
                continue
            if not line: # remove empty lines
                continue 
            words = line.split()
            if words[0].lower() == 'case_num':
                casenum  = int(words[1])
            if words[0].lower() == 'clean_directories':
                if words[1].lower() == 'n':
                    clear_all = 0

    return casenum,clear_all
#------------------------------------------------------------------

# Generate new case number
def generate_new_casenum():
    ldir = glob.glob('casenum_*/')
    if ldir == []:
        casenum = 1
    else:
        casenum = -1
        for dname in ldir:
            spl  = re.split('[_ /]',dname)
            cnum = spl[len(spl)-2]
            if float(cnum) > casenum:
                casenum = float(cnum)
    return int(casenum)+1
#------------------------------------------------------------------

# Add/Change polydispersity values
def editinpfyle(currdir,inpfyle,dispval,runval,nresids,nchains,cnum):

    outfyle = 'pdi_' +str(dispval)+'_run_'+str(runval)+'_'+inpfyle
    fr  = open(inpfyle,'r')
    fw  = open(outfyle,'w')
    fid = fr.read().replace("py_dispinp",str(dispval)).\
              replace("py_resinp",str(nresids)).\
              replace("py_chaininp",str(nchains)).\
              replace("py_casenum",str(cnum))
    fw.write(fid)
    fw.close()
    fr.close()

    return outfyle
#------------------------------------------------------------------

# Run genconf.py from source directory
def run_genconf(inp_fyle,lig_fyles,casenum):
    for fname in lig_fyles:
        if not os.path.exists(fname):
            raise RuntimeError(fname + " not found!")
    print("Running casenumber: ", casenum, "from source")
    subprocess.call(["python", "genconf.py", inp_fyle])
#------------------------------------------------------------------

# Run steps-1,2,3.tcl
def run_all_steps(init_dir):
    os.chdir(init_dir)
    if glob.glob('*.tcl') == []:
        raise RuntimeError("No tcl files not found in "+init_dir)

    subprocess.call(["vmd", "-dispdev", "text", "-e", "step1.tcl"])
    subprocess.call(["vmd", "-dispdev", "text", "-e", "step2.tcl"])
    subprocess.call(["vmd", "-dispdev", "text", "-e", "step3.tcl"])
#------------------------------------------------------------------

# Copy L.pdb/L.psf/*.top to relevant directory
def cpy_pdb_top(init_dir,findir,biomass):
    os.chdir(init_dir)
    # Copy topology
    if glob.glob('*.top') == []:
        raise RuntimeError("No top files formed")
    elif len(glob.glob('*.top')) == 1:
        top_fname = glob.glob(init_dir + '/*.top')
    else:
        fnames = glob.glob(init_dir+'/*.top')
        top_fname = max(fnames, key=os.path.getctime)

    spl_str  = top_fname.split("/")
    src_top  = spl_str[len(spl_str)-1]
    top_file = biomass + "_" + src_top

    # Copy pdb
    if glob.glob('*.pdb') == []:
        raise RuntimeError("No pdb files formed")
    elif len(glob.glob('*.pdb')) == 1:
        pdb_fname = glob.glob(init_dir + '/*.pdb')
    else:
        fnames = glob.glob(init_dir+'/*.pdb')
        pdb_fname = max(fnames, key=os.path.getctime)

    spl_str  = pdb_fname.split("/")
    src_pdb  = spl_str[len(spl_str)-1]
    pdb_file = biomass + "_" + src_pdb

    gencpy2(init_dir,findir,src_top,top_file)
    gencpy2(init_dir,findir,src_pdb,pdb_file)
#------------------------------------------------------------------
# Copy supplementary files
def cpy_supp_files(init_dir,findir,dispval,remflag):
    os.chdir(init_dir)
    # Copy all top files
    list_fnames = glob.glob('*.top')
    if not list_fnames == []:
        for i in range(len(list_fnames)):
            gencpy(init_dir,findir,list_fnames[i])

    # Copy all inp files
    list_fnames = glob.glob('*.inp')
    if not list_fnames == []:
        for i in range(len(list_fnames)):
            gencpy(init_dir,findir,list_fnames[i])

    # Copy all tcl files
    list_fnames = glob.glob('*.tcl')
    if not list_fnames == []:
        for i in range(len(list_fnames)):
            gencpy(init_dir,findir,list_fnames[i])


    # Copy all prm files
    list_fnames = glob.glob('*.prm')
    if not list_fnames == []:
        for i in range(len(list_fnames)):
            gencpy(init_dir,findir,list_fnames[i])


    # Copy all dat files
    list_fnames = glob.glob('*.dat')
    if not list_fnames == []:
        for i in range(len(list_fnames)):
            gencpy(init_dir,findir,list_fnames[i])

    # Copy all log files
    list_fnames = glob.glob('*.log')
    if not list_fnames == []:
        for i in range(len(list_fnames)):
            gencpy(init_dir,findir,list_fnames[i])

    # Copy all log files
    list_fnames = glob.glob('log*')
    if not list_fnames == []:
        for i in range(len(list_fnames)):
            gencpy(init_dir,findir,list_fnames[i])
    
    # Copy all psf files
    list_fnames = glob.glob('*.psf')
    if not list_fnames == []:
        for i in range(len(list_fnames)):
            gencpy(init_dir,findir,list_fnames[i])

    # Copy all pdb files
    list_fnames = glob.glob('*.pdb')
    if not list_fnames == []:
        for i in range(len(list_fnames)):
            gencpy(init_dir,findir,list_fnames[i])

    ## Here on just a contingency at the expense of computational
    # time. Comment when necessary.
    # Copy all txt files
    list_fnames = glob.glob('*.txt')
    if not list_fnames == []:
        for i in range(len(list_fnames)):
            gencpy(init_dir,findir,list_fnames[i])
    
    # Copy all NAMD files
    list_fnames = glob.glob('out*')
    if not list_fnames == []:
        for i in range(len(list_fnames)):
            gencpy(init_dir,findir,list_fnames[i])
            if remflag: #remove for saving space in $HOME
                os.remove(list_fnames[i])

    list_fnames = glob.glob('*.namd')
    if not list_fnames == []:
        for i in range(len(list_fnames)):
            gencpy(init_dir,findir,list_fnames[i])

    list_fnames = glob.glob('*.conf')
    if not list_fnames == []:
        for i in range(len(list_fnames)):
            gencpy(init_dir,findir,list_fnames[i])
#------------------------------------------------------------------
