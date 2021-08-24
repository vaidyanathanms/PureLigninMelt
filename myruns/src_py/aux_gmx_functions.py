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
import datetime
#------------------------------------------------------------------

# General copy script
def gencpy(dum_maindir,dum_destdir,fylname):

    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        raise RuntimeError('ERROR: ' + srcfyl + ' not found')

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
        conf_fname = max(fnames, key=os.path.getmtime)

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
        topol_fname = max(fnames, key=os.path.getmtime)
    
    return conf_fname, topol_fname
#------------------------------------------------------------------

# Count and display number of chains. Assumes dimer as the smallest
# possible polymer.
def count_nchains(destdir,inpname,inp_nchains):

    if not os.path.exists(inpname):
        raise RuntimeError(inpname," not found in ",destdir)
    
    print(inpname)
    ch_id   = 0 # chain ID to distinguish chain beginning
    refmon_num = -1 # dummy value
    mons_perchain_list = []; atoms_perchain_list = []
    fwrlist = open(destdir + '/chainlist.dat','w')
    fwrlist.write('%s\t%s\t%s\n' %('Chain_ID','Num_mons','Num_atoms'))
    with open(inpname) as finp:
        finp.readline() #skip first line
        tot_natoms = int(finp.readline()) # total num atoms
        for line in finp:
            line = line.rstrip('/n')
            if line.startswith('#'): #skip lines with #
                continue
            elif not line: #skip empty lines
                continue
            words=line.split()

            if len(words) == 3: #EOF
                if ch_id == 0:
                    raise RuntimeError('ERROR: Check gro file contents')
                mons_perchain_list.append(mon_cnt)
                atoms_perchain_list.append(atom_cnt)
                fwrlist.write('%d\t%d\t%d\n' %(ch_id,mon_cnt\
                                               ,atom_cnt))
                ch_id += 1
                continue

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
              sum(atoms_perchain_list))
        raise RuntimeError('Error in counting atoms. Please report!')

    num_chains = len(mons_perchain_list)
    if num_chains != inp_nchains:
        print('ERR: num of chains in gro file does not match with'
              'input number of chains', num_chains, inp_nchains) 
        raise RuntimeError('Error - check input number of chains')

    return mons_perchain_list, atoms_perchain_list
#------------------------------------------------------------------

# Retrieve ID
def retrieve_id(inpstr):

    match = re.findall(r'[A-Za-z]+|\d+', inpstr)
    if match == []:
        print('Unknown atom ID: ', inpstr)
        raise RuntimeError('Error - check gro file')
    elif not match[0].isdigit(): 
        print('Unknown atom ID: ', inpstr)
        raise RuntimeError('Error - check gro file')

    return int(match[0])
#------------------------------------------------------------------  

# Find trajector and tpr files
def find_tpr_trr_files(dum_inpdir):
    # check tpr files (*.tpr)
    if glob.glob(dum_inpdir+'/*.tpr') == []:
        print('ERR: No tpr files found in ', dum_inpdir)
        return -1, -1
    elif len(glob.glob(dum_inpdir+'/*.tpr')) == 1:
        conf_fname = glob.glob(dum_inpdir+'/*.tpr')[0]
    elif len(glob.glob(dum_inpdir+'/*.tpr')) > 1:
        print('More than one tpr file found. Using latest')
        fnames = glob.glob(dum_inpdir+'/*.tpr')
        tpr_fname = max(fnames, key=os.path.getmtime)

    # check trr files (*.trr)
    if glob.glob(dum_inpdir+'/*.trr') == []:
        print('ERR: No trr files found in ', dum_inpdir)
        return -1, -1
    elif len(glob.glob(dum_inpdir+'/*.trr')) == 1:
        conf_fname = glob.glob(dum_inpdir+'/*.tpr')[0]
    elif len(glob.glob(dum_inpdir+'/*.trr')) > 1:
        print('More than one trr file found. Using latest')
        fnames = glob.glob(dum_inpdir+'/*.trr')
        trr_fname = max(fnames, key=os.path.getmtime)
        
    return tpr_fname, trr_fname
#------------------------------------------------------------------    
# Create index files for GROMACS
def create_anagrps_inp(destdir,monlist,atomlist,nchains):

    fname = destdir + '/chaininp.inp'
    frglist = open(fname,'w')
    mon_init = 1
    for ncnt in range(nchains):
        mon_fin = mon_init + monlist[ncnt] - 1
        if ncnt != nchains - 1:
            frglist.write("%s %d %s %d;\n" %("resindex", mon_init,\
                                              "to", mon_fin))
        else:
            frglist.write("%s %d %s %d\n" %("resindex", mon_init,\
                                             "to", mon_fin))

        mon_init = mon_fin + 1
    
    frglist.close()        
#------------------------------------------------------------------    
# Run analysis
def run_analysis(nchains,rgflag,msdflag,rdfflag,segrgflag,headdir,\
                 destdir,trajfile,tprfile,conffile,temp):

    if rgflag:
        fname = 'rgcomp_pyinp.sh'
        replace_strings(nchains,headdir,destdir,fname,trajfile,tprfile,conffile,temp,'rg')
    if msdflag:
        fname = 'msdcomp_pyinp.sh'
        replace_strings(nchains,headdir,destdir,fname,trajfile,tprfile,conffile,temp,'msd')
    if rdfflag:
        fname = 'rdfcomp_pyinp.sh'
        replace_strings(nchains,headdir,destdir,fname,trajfile,tprfile,conffile,temp,'rdf')
    if segrgflag:
        fname = 'segment_rgcomp_pyinp.sh'
        replace_strings(nchains,headdir,destdir,fname,trajfile,tprfile,conffile,temp,'segment_rgcomp_pyinp')
#--------------------------------------------------------------------

# Replace strings in job files
def replace_strings(nchains,headdir,destdir,fname,trajfile,tprfile,conffile,temp,anastr):

    if not os.path.exists(headdir + '/' + fname):
        raise RuntimeError(fname + ' not found in ' +headdir)

    jobname  = anastr + '_T_' + str(temp)
    trajfile = split_and_return_filename(trajfile)
    tprfile  = split_and_return_filename(tprfile)
    conffile = split_and_return_filename(conffile)
    
    gencpy(headdir,destdir,fname)
    rev_fname = fname.replace('_pyinp','')
    fr  = open(destdir + '/' + fname,'r')
    fw  = open(destdir + '/' + rev_fname,'w')
    fid = fr.read().replace("py_nchains",str(nchains)).\
          replace("py_trajfile",trajfile).\
          replace("py_tprfile",tprfile).\
          replace("py_conffile",conffile).\
          replace("py_jname",jobname)
    fw.write(fid)
    fw.close()
    fr.close()
#--------------------------------------------------------------------

# Split and return last value from strings for finding files
def split_and_return_filename(inp_fylename):

    outstr = inp_fylename.split('/')
    return outstr[len(outstr)-1]
#--------------------------------------------------------------------

# if __name__ 
if __name__ == 'main':
   main()
#--------------------------------------------------------------------
