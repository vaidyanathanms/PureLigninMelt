# Auxiliary data for generic plot
#------------------------------------------------------------------
# Make labels and graphics
# Version: Mar-15-2021
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
def gencpy(dum_maindir,dum_destdir,fylname):

    srcfyl = dum_maindir + '/' + fylname

    if not os.path.exists(srcfyl):
        print('ERROR: ', srcfyl, 'not found')
        return

    desfyl = dum_destdir + '/' + fylname
    shutil.copy2(srcfyl,desfyl)
#------------------------------------------------------------------

def ret_temp_dir(scr_dir,inp_type,biomass,dispval,casenum,subdir,\
                 solv_type='None'):
    # Check all directories
    head_dir = scr_dir + '/' + inp_type
    if not os.path.isdir(head_dir):
        raise RuntimeError(head_dir + " does not exist!")

       
    poly_dir = head_dir + '/' + biomass
    if not os.path.isdir(poly_dir):
        raise RuntimeError(poly_dir + " does not exist!")
        
    if inp_type == 'melts':
        fig_dir  = poly_dir + '/results_figs'
        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)

        poly_dir = poly_dir + '/pdi_' + str(dispval)
        if not os.path.isdir(poly_dir):
            raise RuntimeError(poly_dir + " does not exist!")

    rundir = poly_dir + '/run_' + str(casenum)
    if not os.path.isdir(rundir):
        raise RuntimeError(rundir + " does not exist!")

    if inp_type == 'solvents':
        workdir1 = rundir + '/' + solv_type
    elif inp_type == 'cosolvents':
        h_workdir = rundir + '/' + solv_type
        workdir1 = h_workdir + '/' + solv_type + '_water'
    else:
        workdir1 = rundir #even for inp_type = melts

    if not os.path.isdir(workdir1):
        raise RuntimeError(workdir1, " does not exist!")

    anadir = workdir1 + '/T_' + str(subdir)
    if not os.path.isdir(workdir1):
        raise RuntimeError(anadir, " does not exist!")

    return workdir1, anadir, fig_dir
#------------------------------------------------------------------        

def ret_ana_dir(scr_dir,inp_type,biomass,dispval,casenum,subdir,\
                 solv_type='None'):
    # Check all directories
    head_dir = scr_dir + '/' + inp_type
    if not os.path.isdir(head_dir):
        raise RuntimeError(head_dir + " does not exist!")

    fig_dir  = head_dir + '/results_figs'
    if not os.path.isdir(fig_dir):
        os.mkdir(fig_dir)
        
    poly_dir = head_dir + '/' + biomass
    if not os.path.isdir(poly_dir):
        raise RuntimeError(poly_dir + " does not exist!")
        
    if inp_type == 'melts':
        poly_dir = poly_dir + '/pdi_' + str(dispval)
        if not os.path.isdir(poly_dir):
            raise RuntimeError(poly_dir + " does not exist!")

    rundir = poly_dir + '/run_' + str(casenum)
    if not os.path.isdir(rundir):
        raise RuntimeError(rundir + " does not exist!")

    if inp_type == 'solvents':
        workdir1 = rundir + '/' + solv_type
    elif inp_type == 'cosolvents':
        h_workdir = rundir + '/' + solv_type
        workdir1 = h_workdir + '/' + solv_type + '_water'
    else:
        workdir1 = rundir #even for inp_type = melts

    if not os.path.isdir(workdir1):
        raise RuntimeError(workdir1, " does not exist!")

    anadir = workdir1 + '/' + subdir
    if not os.path.isdir(workdir1):
        raise RuntimeError(anadir, " does not exist!")

    return workdir1, anadir, fig_dir
#------------------------------------------------------------------        
# Change width of seaborn barplot
def change_width(ax, new_value):
    for patch in ax.patches:
        current_width = patch.get_width()
        diff = current_width - new_value

        # Change the bar width
        patch.set_width(new_value)

        # Recenter the bar
        patch.set_x(patch.get_x() + diff * .5)
#------------------------------------------------------------------

# Return Chain ID from file name
def ret_chid(inpfyle):
    outstr = re.split(r'/|.xvg|_',inpfyle)
    chid = outstr[len(outstr)-2]
    if not chid.isdigit():
        print('ERR: not a numeric value', chid)
    else:
        return int(chid)
#------------------------------------------------------------------

def ret_mons(inpfyle):
    if not os.path.exists(inpfyle):
        raise RuntimeError('File not found ' + inpfyle)
    mon_arr = []
    with open(inpfyle) as finp:
        finp.readline()
        for line in finp:
            line = line.rstrip('/n')
            outstr = line.split()
            mon_arr.append(outstr[len(outstr)-2])
    return mon_arr
#------------------------------------------------------------------
# All Rg plots
def plot_allrg(df,fig_dir,pdi_val):

    #Plot Rg^2 - T
    maxval = df['<Rg^2>'].max()
    figa, axa = plt.subplots()
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()
    axa.scatter(x=df['Temperature'],y=df['<Rg^2>'])
    axa.set_ylabel(r'Temperature ($K$)')
    axa.set_ylabel(r'$R_{g}^2$ ($nm^2$)')
    change_width(axa,0.2)
    figa.savefig(fig_dir + '/'+'rg2_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)
    
    #Plot Rg^4 - T
    maxval = df['<Rg^4>'].max()
    figa, axa = plt.subplots()
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()
    axa.scatter(x=df['Temperature'],y=df['<Rg^4>'])
    axa.set_ylabel(r'Temperature ($K$)')
    axa.set_ylabel(r'$R_{g}^4$ ($nm^4$)')
    change_width(axa,0.2)
    figa.savefig(fig_dir + '/'+'rg4_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)
    
    #Plot Alpha - T
    maxval = df['Rg4/Rg2^2'].max()
    figa, axa = plt.subplots()
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()
    axa.scatter(x=df['Temperature'],y=df['Rg4/Rg2^2'])
    axa.set_ylabel(r'Temperature ($K$)')
    axa.set_ylabel(r'$<Rg^4>/<Rg^2>^2$')
    change_width(axa,0.2)
    figa.savefig(fig_dir + '/'+'alpha' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)
#------------------------------------------------------------------
# Plot Tg
def plot_tg(df,fig_dir,pdi_val):
    maxval = df['SV_NPT'].max()
    figa, axa = plt.subplots()
    plt.style.use('seaborn-colorblind')
    plt.tight_layout()
    axa.scatter(x=df['Temp'],y=df['SV_NPT'])
    axa.set_xlabel(r'Temperature ($K$)')
    axa.set_ylabel(r'Specific Volume ($cm^3/g$)')
    change_width(axa,0.2)
    figa.savefig(fig_dir + '/'+'sv_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    return figa, axa
#------------------------------------------------------------------

# if __name__
if __name__ == '__main__':
    main()
#------------------------------------------------------------------
