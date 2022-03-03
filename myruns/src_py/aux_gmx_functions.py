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
import parmed as pmd
import general_functions as genfun
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
    # check structure files (.gro)
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
        print("WARNING: No polymer topology files found")
        topol_fname = 'None'
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
        tpr_fname = glob.glob(dum_inpdir+'/*.tpr')[0]
    elif len(glob.glob(dum_inpdir+'/*.tpr')) > 1:
        print('More than one tpr file found. Using latest')
        fnames = glob.glob(dum_inpdir+'/*.tpr')
        tpr_fname = max(fnames, key=os.path.getmtime)

    # check trr files (*.trr)
    if glob.glob(dum_inpdir+'/*.trr') == []:
        print('ERR: No trr files found in ', dum_inpdir)
        return -1, -1
    elif len(glob.glob(dum_inpdir+'/*.trr')) == 1:
        trr_fname = glob.glob(dum_inpdir+'/*.tpr')[0]
    elif len(glob.glob(dum_inpdir+'/*.trr')) > 1:
        print('More than one trr file found. Using latest')
        fnames = glob.glob(dum_inpdir+'/*.trr')
        trr_fname = max(fnames, key=os.path.getmtime)
        
    return tpr_fname, trr_fname
#------------------------------------------------------------------    
# Create index files for GROMACS: Rg/ShapeFactor/SegRg
def create_rggrps_inp(destdir,monlist,atomlist,nchains):

    fname = destdir + '/rgchaininp.inp'
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
# Create index files for GROMACS: RDFs
def create_rdfgrps_inp(destdir,monlist,atomlist,nchains):

    #Inter/Intrachain RDF
    mon_init = 1
    for ncnt in range(nchains):
        fname = destdir + '/rdf' + str[ncnt+1] + '.inp'
        frdflist = open(fname,'w')
        mon_fin = mon_init + monlist[ncnt] - 1
        frdflist.write("%s %d %s %d;\n" %("name C1 and resindex",\
                                          mon_init,"to", mon_fin))
        frdflist.write("%s %d %s %d;\n" %("name C1 and resindex",\
                                          mon_init,"to", mon_fin))
        frdflist.write("%s %d %s %d\n" %("name C1 and not resindex"\
                                          ,mon_init,"to", mon_fin))
        frdflist.close()        
        mon_init = mon_fin + 1

    #HH/HG/HS
    with open(destdir + '/Hall.inp') as frdflist:
        frdflist.write("%s;\n" %("resname PHP and name C1"))
        frdflist.write("%s;\n" %("resname PHP and name C1"))
        frdflist.write("%s;\n" %("resname GUAI and name C1"))
        frdflist.write("%s\n" %("resname SYR and name C1"))
        
    #GH/GG/GS
    with open(destdir + '/Hall.inp') as frdflist:
        frdflist.write("%s;\n" %("resname GUAI and name C1"))
        frdflist.write("%s;\n" %("resname PHP and name C1"))
        frdflist.write("%s;\n" %("resname GUAI and name C1"))
        frdflist.write("%s\n" %("resname SYR and name C1"))

    #SH/SG/SS
    with open(destdir + '/Hall.inp') as frdflist:
        frdflist.write("%s;\n" %("resname SYR and name C1"))
        frdflist.write("%s;\n" %("resname PHP and name C1"))
        frdflist.write("%s;\n" %("resname GUAI and name C1"))
        frdflist.write("%s\n" %("resname SYR and name C1"))

    #FF/FP/FG/FS
    with open(destdir + '/Hall.inp') as frdflist:
        frdflist.write("%s;\n" %("resname FERUT and name C1"))
        frdflist.write("%s;\n" %("resname FERUT and name C1"))
        frdflist.write("%s;\n" %("resname PCA and name C1"))
        frdflist.write("%s;\n" %("resname GUAI and name C1"))
        frdflist.write("%s\n" %("resname SYR and name C1"))

    #PF/PP/PG/PS
    with open(destdir + '/Hall.inp') as frdflist:
        frdflist.write("%s;\n" %("resname PCA and name C1"))
        frdflist.write("%s;\n" %("resname FERUT and name C1"))
        frdflist.write("%s;\n" %("resname PCA and name C1"))
        frdflist.write("%s;\n" %("resname GUAI and name C1"))
        frdflist.write("%s\n" %("resname SYR and name C1"))
#------------------------------------------------------------------ 

# Create index files for GROMACS: HBs
def create_hbgrps_inp(destdir,monlist,atomlist,nchains):
    #Inter/Intrachain HB
    mon_init = 1
    for ncnt in range(nchains):
        fname = destdir + '/hb' + str[ncnt+1] + '.inp'
        fhblist = open(fname,'w')
        mon_fin = mon_init + monlist[ncnt] - 1
        fhblist.write("%s %d %s %d;\n" %("resindex",mon_init,\
                                          "to", mon_fin))
        fhblist.write("%s %d %s %d;\n" %("resindex",mon_init,\
                                          "to", mon_fin))
        fhblist.write("%s %d %s %d\n" %("not resindex",mon_init,\
                                        "to", mon_fin))
        frdflist.close()        
        mon_init = mon_fin + 1

    #FA-FA/FA-PCA/FA-G/FA-S
    with open(destdir + '/resinp.inp') as flist:
        flist.write("%s;\n" %("resname FERUT"))
        flist.write("%s;\n" %("resname PCA"))
        flist.write("%s;\n" %("resname PHP"))
        flist.write("%s;\n" %("resname GUAI"))
        flist.write("%s\n" %("resname SYR"))
#------------------------------------------------------------------ 

# Run analysis - OBSOLETE on Feb-14-2022. Added to gmx_functions.py
def run_analysis(nchains,rgflag,msdflag,rdfflag,segrgflag,shapeflag\
                 ,headdir,destdir,trajfile,tprfile,conffile,temp):

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
    if shapeflag:
        fname = 'segment_rgcomp_pyinp.sh'
        replace_strings(nchains,headdir,destdir,fname,trajfile,tprfile,conffile,temp,'shapecomp_pyinp')
#--------------------------------------------------------------------

# Replace strings and set job files
def set_jobfile(nchains,headdir,destdir,fname,trajfile,tprfile,conffile,temp,anastr):

    if not os.path.exists(headdir + '/' + fname):
        raise RuntimeError(fname + ' not found in ' +headdir)
  
    gencpy(headdir,destdir,fname)
    rev_fname = fname.replace('_pyinp','')
    fr  = open(destdir + '/' + fname,'r')
    fw  = open(destdir + '/' + rev_fname,'w')
    if anastr == 'rho':
        jobname  = anastr
        fid = fr.read().replace("py_nchains",str(nchains)).\
              replace("py_jname",jobname).\
              replace("py_Tinit",str(temp[0])).\
              replace("py_Tfin",str(temp[1])).\
              replace("py_dT",str(temp[2]))
    else:
        jobname  = anastr + '_T_' + str(temp)
        trajfile = split_and_return_filename(trajfile)
        tprfile  = split_and_return_filename(tprfile)
        conffile = split_and_return_filename(conffile)
        fid = fr.read().replace("py_nchains",str(nchains)).\
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

# Check for PSF files
def check_psf(inpdir):
    desdir = inpdir + '/init_files'
    if not os.path.isdir(desdir):
        print(desdir + " does not exist!")
        os.mkdir(desdir)
        return -1
    else:
        # Default PSF file is L.psf. Just to be consistent I recreate
        # PSF file using ParmEd - compRg.psf such that the resnames are
        # FERU and PCA. Also, while calculating the Rg, the "fragment" 
        # identifier is used as opposed to segname identifier.
        if not os.path.exists(desdir+'/compRg.psf'):
            print('compRg.psf for segmental Rg not found in ' +  desdir)
            print('Will generate compRg.psf files using top/gro files')
            return -1
        else:
            return 1
#--------------------------------------------------------------------

# Copy relevant gro/top files
def cpy_gro_top_files(inpdir,desdir):
    grofile = genfun.extract_file_name(genfun.find_latest_file(inpdir,'gro'))
    topfile = genfun.extract_file_name(genfun.find_latest_file(inpdir,'top'))
    gencpy(inpdir,desdir,grofile)
    gencpy(inpdir,desdir,topfile)
    convert_top_to_psf(topfile,grofile,desdir)
#--------------------------------------------------------------------

# Convert top to psf/crd
# Thank you, Alan Hicks :) for introducing to ParmEd
# https://github.com/ParmEd/ParmEd
def convert_top_to_psf(inptop,inpgro,destdir):
    print("Generating compRg.psf file...")
    gmx_top = pmd.load_file(destdir + '/' + inptop, \
                            xyz= destdir + '/' + inpgro)

    # ParmEd doesn't like overwriting. So delete and recreate
    if os.path.exists(destdir + '/inp_amber.top'):
        os.remove(destdir + '/inp_amber.top')
    if os.path.exists(destdir + '/inp_amber.crd'):
        os.remove(destdir + '/inp_amber.crd')

    gmx_top.save(destdir + '/inp_amber.top', format='amber')
    gmx_top.save(destdir + '/inp_amber.crd', format='rst7')   
    amber = pmd.load_file(destdir + '/inp_amber.top',\
                          destdir + '/inp_amber.crd')
    amber.save(destdir + '/compRg.psf') #for computing segment Rg
    amber.save(destdir + '/compRg.crd')
#--------------------------------------------------------------------

# if __name__ 
if __name__ == 'main':
   main()
#--------------------------------------------------------------------
