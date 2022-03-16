# Analyze and plot density/Tg
# Version: Mar-16-2022
#------------------------------------------------------------------

# Import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
import sys
import re
import shutil
import glob
import math
import subprocess
from plt_aux import *
from compute_props import *
#------------------------------------------------------------------
# Input data for analysis
pdi_val   = 1.8
run_arr   = [1] # run numbers for a given biomass
temp_min  = 250 # Minimum temperature
temp_max  = 501 # Maximum temperature
temp_dt   = 10  # Temperature dt
pdi_arr   = [1.0,1.8,3.0]
mark_arr  = ['o','d','s']
nchains   = 20
#------------------------------------------------------------------
# Global arrays
denarr = np.arange(temp_min,temp_max,5*temp_dt) #plotting time data
#------------------------------------------------------------------
# Plot glass transition for all systems
print("Analyzing Tg data")
for bio_indx in range(len(biom_arr)): # loop in biomass
    
    biomass = biom_arr[bio_indx]
    res_dir = scr_dir + '/' + inp_type + '/'+ biomass +\
              '/pdi_' + str(pdi_val) + '/results'
    if not os.path.isdir(res_dir):
        os.mkdir(res_dir)
        
    # Write data and then read to concatenate
    tg_fyl = res_dir + '/tgdata.dat'
    fc_tg = open(tg_fyl,'w')
    fc_tg.write('%s\t%s\t%s\n' %('Biomass','Temp','SV_NPT'))
    
    tg_avg = np.zeros(0)

    # Define axes labels for distribution
    fig1,ax1 = plt.subplots()
    set_axes(ax1,plt,r'Time (ps)',r'Density ($kg$/$m^3$)')
    
    yrho_avg = np.zeros(0) # Average density

    for tval in range(temp_min,temp_max,temp_dt): # loop in temp
        temp_leg  = str(tval)
        yall = np.zeros(0)
        yavg = 0
        
        for casenum in range(len(run_arr)): # loop in runarr
            wdir,tdir,fig_dir = ret_temp_dir(scr_dir,inp_type,biomass,\
                                             pdi_val,run_arr[casenum],\
                                             tval,solv_type)
            print("Analyzing: ", pdi_val,tval,run_arr[casenum])
            # Check file
            if os.path.exists(tdir + '/dens_npt.xvg'):
                fname  = tdir + '/dens_npt.xvg'
            elif os.path.exists(wdir + '/dens_npt.xvg'):
                fname  = wdir + '/dens_npt.xvg'
            else:
                print(fname, " does not exist! ")
                print("ERR: Compute densities before Tg calculation..")
                continue
                
            # Open and parse file
            with open(fname) as fin:
                lines = (line.lstrip() for line in fin \
                         if not line.lstrip().startswith('#') and \
                         not line.lstrip().startswith('@'))
                data  = np.loadtxt(lines)

            yall = np.append(yall,data[:,1]) #append all y-data
            startpt = int(0.1*len(yall))
            yavg += np.average(data[startpt:len(yall)-1,1])
            
            if casenum == 0 and tval in denarr: # plot time-density variations
                ax1.plot(data[:,0],data[:,1],label=temp_leg)
            
        # append avg to yrho_avg
        yrho_avg = np.append(ytemp_avg,yavg/len(run_arr))
        
    # Save density-time plots
    ax1.legend(loc=0)
    fig1.savefig(fig_dir + '/'+'denstime_'+str(pdi_val)+'.png',\
                 dpi=fig1.dpi)
    plt.close(fig1)

    # Plot all data
    indx = 0
    for t_plt in range(temp_min,temp_max,temp_dt): # loop in temp
        if yrho_avg[indx] == 0:
            continue
            indx += 1
        sval = 1000.0/yrho_avg[indx]
        fc_tg.write('%s\t%g\t%g\n' %(biomass,t_plt,sval))
        indx += 1

    fc_tg.close()
    df=pd.read_table(tg_fyl)
    figa, axa = plot_tg(df,fig_dir,pdi_val) 
    tgval = compute_tg(df,axa) #compute tgval
    figa.savefig(fig_dir + '/'+'sv2_' + str(pdi_val) + '.png',\
                 dpi=figa.dpi)
    plt.close(figa)
    
    if len(pdi_arr) > 1: #plot as a function of PDI

        fig2, ax2 = plt.subplots()
        set_axes(ax2,plt,r'Temperature ($K$)',r'Specific Volume ($cm^{3}/g$)')

        for pdiindx in range(len(pdi_arr)):
            pdiplt = pdi_arr[pdiindx]
            pdileg = 'PDI: ' + str(pdiplt)
            res_dir = scr_dir + '/' + inp_type + '/'+ biomass +\
                      '/pdi_' + str(pdiplt) + '/results'

            if not os.path.isdir(res_dir):
                print(res_dir, "does not exist")
                continue

            else:
                tg_fyl = res_dir + '/tgdata.dat'
                df=pd.read_table(tg_fyl)
                print('Plotting', pdileg)
                ax2.scatter(x=df['Temp'],y=df['SV_NPT'],label=pdileg)

        fig2.savefig(fig_dir + '/'+'sv_all.png',dpi=fig2.dpi)
        plt.close(fig2)
#------------------------------------------------------------------
