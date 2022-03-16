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
from plt_inps import *
from plt_aux import *
from compute_props import *
#------------------------------------------------------------------

# Input data
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
# Plot RDF (pol-solvent, pol-water)
if rdf_pol:
    
    print("Analyzing whole polymer RDF data")
    
    for bio_indx in range(len(biom_arr)): # loop in biomass
        biomass = biom_arr[bio_indx]
        yrg_avg = np.zeros(0)

        # Plot polymer-solvent RDF
        fig1,ax1 = plt.subplots()
        ax1.set_xlabel(r'$r$ (nm)')
        ax1.set_ylabel(r'$g_{\mathrm{pol-orgsol}}(r)$')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()

        # Plot polymer-water RDF
        fig2,ax2 = plt.subplots()
        ax2.set_xlabel(r'$r$ (nm)')
        ax2.set_ylabel(r'$g_{\mathrm{pol-water}}(r)$')
        plt.style.use('seaborn-colorblind')
        plt.tight_layout()

        for sol_indx in range(len(otyp_arr)): # loop in systems
            solv_type = otyp_arr[sol_indx]
            solv_leg  = otyp_leg[sol_indx]
            xdata = np.zeros(0)
            y1data = np.zeros(0); y2data = np.zeros(0)
            rdata = np.zeros(0) # for counting repeats

            for casenum in range(len(run_arr)): # loop in runarr
                wdir,anadir,fig_dir = ret_ana_dir(scr_dir,inp_type,biomass,\
                                                  disperse,run_arr[casenum],\
                                                  'anafiles',solv_type)
                print("Analyzing: ", biomass,solv_type,run_arr[casenum])
                # Check file
                if os.path.exists(anadir + '/rdfout1.xvg'):
                    fname  = anadir + '/rdfout1.xvg'
                elif os.path.exists(wdir + '/rdfout1.xvg'):
                    fname  = wdir + '/rdfout1.xvg'
                else:
                    print("rdfout1.xvg does not exist! ")
                    continue

                # Open and parse file
                with open(fname) as fin:
                    lines = (line.lstrip() for line in fin \
                             if not line.lstrip().startswith('#') and \
                             not line.lstrip().startswith('@'))
                    data  = np.loadtxt(lines)

                if casenum == 0: #append all x/y-data
                    xdata = np.append(xdata,data[:,0])
                    y1data = np.append(y1data,data[:,1])
                    y2data = np.append(y2data,data[:,2])
                    omax_xval = np.amax(xdata)
                    olen_data = len(xdata)
                    rdata = np.full((len(xdata)),1,dtype=int)
                else:
                    if np.amax(data[:,0]) >= omax_xval:
                        for i_indx in range(olen_data):
                            y1data[i_indx] += data[i_indx,1]
                            y2data[i_indx] += data[i_indx,2]
                            rdata[i_indx] += 1
                        len_data = len(data[:,0])
                        max_xval = np.amax(data[:,0])
                        xdata = np.append(xdata,data[i_indx:len_data-1,0])
                        y1data = np.append(y1data,data[i_indx:len_data-1,1])
                        y2data = np.append(y2data,data[i_indx:len_data-1,1])
                        rnew  = np.full((1,len_data-olen_data),1,dtype=int)
                        rdata = np.append(rdata,rnew)
                        omax_xval = max_xval
                        olen_data = len_data
                    else:
                        for i_indx in range(len(data[:,0])):
                            y1data[i_indx] += data[i_indx,1]
                            y2data[i_indx] += data[i_indx,2]
                            rdata[i_indx] += 1
                        
            # Divide ydata by rdata
            y1data /= rdata
            y2data /= rdata
            ax1.plot(xdata,y1data,label=solv_leg)
            ax2.plot(xdata,y2data,label=solv_leg)

        ax1.legend(loc=0)
        ax2.legend(loc=0)
        ax1.set_xlim([0,3.3])
        ax2.set_xlim([0,3.3])
        fig1.savefig(fig_dir + '/'+biomass+'_rdfpsol.png',dpi=fig1.dpi)
        fig2.savefig(fig_dir + '/'+biomass+'_rdfpwat.png',dpi=fig2.dpi)
        plt.close(fig1)
        plt.close(fig2)
#------------------------------------------------------------------
