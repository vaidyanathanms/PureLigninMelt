# Analyze and plot density/Tg
# Version: Mar-16-2022
#------------------------------------------------------------------

# Import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
import math
from plt_aux import *
from compute_props import *
from plt_inps import *
#------------------------------------------------------------------

# Input data for analysis
run_arr   = [1,2,3,4] # run numbers for a given biomass
temp_min  = 250 # Minimum temperature
temp_max  = 501 # Maximum temperature
temp_dt   = 10  # Temperature dt
pdi_arr   = [1.0,1.8,3.0,'expts']
mark_arr  = ['o','d','s']
nchains   = 20
anadir_head = 'all_dens' # change this for different analysis
start_frac  = 0.3 # starting point for averaging
#------------------------------------------------------------------

# Global arrays
denarr = np.arange(temp_min,temp_max,5*temp_dt) #plotting time data
#------------------------------------------------------------------

# Plot glass transition for all systems
print("Analyzing Tg data")
anaout_dir = anaout_dir + '/dens_results' # result_outputs
if not os.path.isdir(anaout_dir):
    os.mkdir(anaout_dir)

for pdi_val in pdi_arr:

    # Set/Check directories
    simout_dir = all_dir + '/results_pdi_' + str(pdi_val)
    if not os.path.isdir(simout_dir):
        print("ERR: " +  simout_dir + " does not exist")       
        continue
    
    # Define axes labels for distribution
    fig1,ax1 = plt.subplots()
    set_axes(ax1,plt,r'Time (ps)',r'Density ($kg$/$m^3$)')

    # All outputs together
    fall_out = open(anaout_dir +'/alltgdata_'+str(pdi_val)+'.dat','w')
    fall_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %('Temp','Dens_R1','Dens_R2',\
                                                'Dens_R3','Dens_R4','TotCase'))
    
    # Set arrays
    yrho_avg = [] # Average density
    ncasetot = [] # Total cases per temperature

    # Temperature loops and averaging
    for tval in range(temp_min,temp_max,temp_dt): # loop in temp
        temp_leg  = str(tval)
        yall = []
        yavg = 0; ncases_pertemp = 0
        fall_out.write('%g\t' %(tval))
        
        for casenum in range(len(run_arr)): # loop in runarr
            
            wdir = simout_dir + '/run_' + str(run_arr[casenum]) +\
                '/T_' + str(tval) + '/' + anadir_head
            if not os.path.isdir(wdir):
                fall_out.write('%s\t' %('N/A'))
                print("ERR: " + wdir + " does not exist")
                continue
            
            print("Analyzing: ", pdi_val,tval,run_arr[casenum])
            if os.path.exists(wdir + '/dens_npt.xvg'): #check file exists
                fname  = wdir + '/dens_npt.xvg' 
            else:
                fall_out.write('%s\t' %('N/A'))
                print(fname, " does not exist! ")
                print("ERR: Compute densities before Tg calculation..")
                continue
                
            # Open and parse file
            with open(fname) as fin:
                lines = (line.lstrip() for line in fin \
                         if not line.lstrip().startswith('#') and \
                         not line.lstrip().startswith('@'))
                data  = np.loadtxt(lines)

            yall = np.append(yall,data[:,1]) #append new y-data
            startpt = int(start_frac*len(yall))
            rho_case = np.average(data[startpt:len(yall)-1,1])
            yavg += rho_case
            ncases_pertemp += 1
            fall_out.write('%g\t' %(rho_case)) # Write to file
             
            # plot time-density variations
            if casenum == 0 and tval in denarr: #end of casenum loop
                ax1.plot(data[:,0],data[:,1],label=temp_leg)
            
        # Do NOT continue if zero cases are found
        if ncases_pertemp == 0:
            fall_out.write('%g\n' %(ncases_pertemp))
            continue
        
        # append avg/totcases to yrho_avg
        yrho_avg = np.append(yrho_avg,yavg/ncases_pertemp) 
        ncasetot = np.append(ncasetot,ncases_pertemp) 
        fall_out.write('%g\n' %(ncases_pertemp)) #end of tval loop
        
    # Do NOT continue if no PDI-temps are found
    if len(yrho_avg) == 0: #to account for boolean values
        continue

    fall_out.close() #Close temp file
    
    # Save density-time plots
    ax1.legend(loc=0)
    fig1.savefig(figout_dir+'/'+'denstime_'+str(pdi_val)+'.png',dpi=fig1.dpi)
    fig1.savefig(figout_dir+'/'+'denstime_'+str(pdi_val)+'.eps',format='eps')
    plt.close(fig1)

    # Write SV-temperature data
    with open(anaout_dir +'/tgdata_'+str(pdi_val)+'.dat','w') as fc_tg:
        fc_tg.write('%s\t%s\t%s\n' %('Temp','N_Cases','SV_NPT'))
        indx = 0
        for t_plt in range(temp_min,temp_max,temp_dt): # loop in temp
            if yrho_avg[indx] == 0:
                indx += 1
                continue
            # Specific Volume: 1000/yrho_avg[indx] (in cm^3/g)
            fc_tg.write('%s\t%g\t%g\n' %(t_plt,ncasetot[indx]\
                                         ,1000.0/yrho_avg[indx]))
            indx += 1

    # Plot SV-temperature data for each PDI
    df=pd.read_table(anaout_dir +'/tgdata_'+str(pdi_val)+'.dat')
    figa,axa = plt.subplots()
    set_axes(axa,plt,r'Temperature ($K$)',r'Specific Volume ($cm^3$/$g$)')
    tgval = compute_tg(df,axa) #compute tgval
    figa.savefig(figout_dir+'/'+'svt_'+str(pdi_val)+'.png',dpi=figa.dpi)
    figa.savefig(figout_dir+'/'+'svt_'+str(pdi_val)+'.eps',format='eps')
    plt.close(figa) # End of PDI loop
    
#------- Plot SV-Temp data for all PDI values together-----------------
fig2, ax2 = plt.subplots()
set_axes(ax2,plt,r'Temperature ($K$)',r'Specific Volume ($cm^{3}/g$)')

for pdi_val in pdi_arr:
    if pdi_val == 'expts':
        pdileg = 'PDI: Experimental Distribution'
    else:
        pdileg = 'PDI: ' + str(pdi_val)
    fname = '/tgdata_'+str(pdi_val)+'.dat'
    if not os.path.exists(anaout_dir + '/' + fname):
        print('ERR: '+fname+' does not exist in ' + anaout_dir)
        continue
    
    df=pd.read_table(anaout_dir + '/' + fname)
    print('Plotting', pdi_val)
    ax2.scatter(x=df['Temp'],y=df['SV_NPT'],label=pdileg)
    
fig2.savefig(figout_dir + '/'+'svt_allpdi.png',dpi=fig2.dpi)
fig2.savefig(figout_dir + '/'+'svt_allpdi.eps',format='eps')
plt.close(fig2)
#------------------------------------------------------------------
