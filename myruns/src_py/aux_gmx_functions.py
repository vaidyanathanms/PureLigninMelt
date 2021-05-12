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
        workdir1 = h_workdir + '/' + solv_type + '_water'
    else:
        workdir1 = rundir #even for inp_type = melts

    if not os.path.isdir(workdir1):
        raise RuntimeError(workdir1, "does not exist")
    
    return workdir1
#------------------------------------------------------------------

#Check pdb/psf/top files for the melt (or polymer)
def check_inp_files(dum_inpdir,top_name):
    # check structure files (*.pdb/.gro)
    if glob.glob(dum_inpdir+'/*.gro') == []:
        raise RuntimeError("No polymer gro files found")
    elif len(glob.glob(dum_inpdir+'/*.gro')) == 1:
        conf_fname = glob.glob(dum_inpdir+'/*.gro')[0]
    elif len(glob.glob(dum_inpdir+'/*.gro')) > 1:
        print('More than one config file found. Using latest')
        fnames = glob.glob(dum_inpdir+'/*.gro')
        conf_fname = max(fnames, key=os.path.getctime)

    # check topology files
    if glob.glob(dum_inpdir+'/*.top') == []:
        raise RuntimeError("No polymer topology files found")
    if top_name != 'None':
        if not os.path.exists(dum_inpdir + '/' + top_name):
            raise RuntimeError("Specified poly top file not found")
        else:
            topol_fname = top_name
    elif len(glob.glob(dum_inpdir+'/*.top')) == 1 and \
         top_name == 'None':
        topol_fname = glob.glob(dum_inpdir+'/*.top')[0]
    else:
        print('More than one topology file found. Using latest')
        fnames = glob.glob(dum_inpdir+'/*.top')
        topol_fname = max(fnames, key=os.path.getctime)
    
    return conf_fname, topol_fname
#------------------------------------------------------------------

# Count and display number of chains. Assumes dimer as the smallest
# possible polymer.
def count_nchains(destdir,inpname,inp_nchains):

    if not os.path.exists(destdir + '/' + inpname):
        raise RuntimeError(inpname," not found in ",destdir)
    
    ch_id   = 0 # chain ID to distinguish chain beginning
    refmon_num = -1 # dummy value
    mons_perchain_list = []; atoms_perchain_list = []
    fwrlist = open(destdir + '/chainlist.dat','w')
    with open(destdir + '/' + inpname) as finp:
        next(finp) # skip first line
        tot_natoms = int(finp.readline()) # total num atoms
        for line in finp:
            line = line.rstrip('/n')
            if line.startswith('#'): #skip lines with #
                continue
            if not line: #skip empty lines
                continue
            words=line.split()
            mon_num = retrieve_id(words[0])

            if mon_num == refmon_num:
                atom_cnt += 1
                continue # skip atoms of the same monomer
            elif mon_num == 1:
                if ch_id != 0: # 0 to initialize
                    mons_perchain_list.append(mon_cnt)
                    atoms_perchain_list.append(atom_cnt)
                    fwrlist.write('%d\t%d\t%d\n' %(ch_id,mon_cnt\
                                                   ,atom_cnt))
                ch_id += 1
                mon_cnt = 1
                atom_cnt = 1
                refmon_num = mon_num
            else:
                mon_cnt += 1
                atom_cnt += 1
                refmon_num = mon_num


    fwrlist.close()                
    # Sanity checks   
    if tot_natoms != sum(atoms_perchain_list):
        print('Total number of atoms do not match', tot_natoms, \
              atoms_perchain_list)
        raise RuntimeError('Error in counting atoms. Please report!')

    num_chains = len(mons_perchain_list)
    if num_chains != inp_nchains:
        print('ERR: num of chains in gro file does not match  with 
        input number of chains', num_chains, inp_nchains) 
        raise RuntimeError('Error - check input number of chains')

    return mons_perchain_list, atoms_perchain_list
#------------------------------------------------------------------

# Retrieve ID
def retrieve_id(inpstr):

    match = re.match(r"([a-z]+)([0-9]+)", inpstr, re.I)
    if match:
        items = match.groups()

    if items == [] or match == []:
        print('Unknown atom ID: ', inpstr)
        raise RunTimeError('Error - check gro file')
    elif not match.isdigit(items[0]): 
        print('Unknown atom ID: ', inpstr)
        raise RunTimeError('Error - check gro file')

    return int(items[0])
#------------------------------------------------------------------    

# Create index files for GROMACS
def create_anaindx_grps(inpfile,monlist,atomlist):


