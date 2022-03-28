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
temp_arr  = [300,400,500]
pdi_arr   = [1.0,1.8,3.0]
mark_arr  = ['o','d','s']
nchains   = 20
anadir_head = 'all_rdfs' # change this for different analysis
#------------------------------------------------------------------
# Plot RDF (intra-RDF, inter-RDF of C1 carbons)
print("Analyzing inter/intra RDF")
anaout_dir = anaout_dir + '/rdf_results' # result_outputs
if not os.path.isdir(anaout_dir):
    os.mkdir(anaout_dir)
    
for pdi_val in pdi_arr:

    # Set/Check directories
    simout_dir = all_dir + '/results_pdi_' + str(pdi_val)
    if not os.path.isdir(simout_dir):
        print("ERR: " +  simout_dir + " does not exist")       

    # Define axes labels for RDFs
    fig1,ax1 = plt.subplots()
    set_axes(ax1,plt,r'$r$ (nm)',r'$g_{\rm{inter}}(r)$')
    fig2,ax2 = plt.subplots()
    set_axes(ax1,plt,r'$r$ (nm)',r'$g_{\rm{intra}}(r)$')
    
    for tval in temp_arr: # loop in temperature
        temp_leg  = str(tval)
        xdata  = np.zeros(0)
        y1data = np.zeros(0) # inter rdf
        y2data = np.zeros(0) # intra rdf
        rdata  = np.zeros(0) # for counting repeats
        
        for casenum in range(len(run_arr)): # loop in runarr
            wdir = simout_dir + '/run_' + str(run_arr[casenum]) +\
                '/T_' + str(tval) + '/' + anadir_head

            if not os.path.isdir(wdir):
                print("ERR: " + wdir + " does not exist")
                continue
            
            print("Analyzing: ", pdi_val,tval,run_arr[casenum])
            flist = glob.glob(wdir + '/rdf_nptmain_*.xvg')
            if not len(flist): #check file(s) exists
                print("ERR: No RDF files exist! ")
                continue
            if len(flist) != nchains:
                print("ERR: Mismatch in number of RDF files ",\
                      len(flist), nchains)
            
            for fname in flist:
                # Open and parse file
                with open(fname) as fin:
                    lines = (line.lstrip() for line in fin \
                             if not line.lstrip().startswith('#') and \
                             not line.lstrip().startswith('@'))
                    data  = np.loadtxt(lines)
                
                if casenum == 0: #append all x/y-data
                    xdata  = np.append(xdata,data[:,0])
                    y1data = np.append(y1data,data[:,1]) #inter
                    y2data = np.append(y2data,data[:,2]) #intra
                    omax_xval = np.amax(xdata)
                    olen_data = len(xdata)
                    rdata = np.full((len(xdata)),1,dtype=int) #bin counts
                else: # to account for different box sizes for each case
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
                            rdata[i_indx] += 1 #end if np.amax >= omax_xval
                # end if casenum == 0
            # end for fname in flist
        #end for casenum in range(len(run_arr))
        
    # Divide ydata by rdata*nchains
    y1data /= (rdata*nchains)
    y2data /= (rdata*nchains) #end for case
    
    # Plot and save RDF plots
    ax1.plot(xdata,y1data,label=solv_leg)
    ax2.plot(xdata,y2data,label=solv_leg)
    
    ax1.legend(loc=0)    
    fig1.savefig(figout_dir+'/'+'interrdf_'+str(pdi_val)+'.png',dpi=fig1.dpi)
    fig1.savefig(figout_dir+'/'+'interrdf_'+str(pdi_val)+'.eps',format='eps')
    plt.close(fig1)

    ax2.legend(loc=0)
    fig2.savefig(figout_dir+'/'+'intrardf_'+str(pdi_val)+'.png',dpi=fig1.dpi)
    fig2.savefig(figout_dir+'/'+'intrardf_'+str(pdi_val)+'.eps',format='eps')
    plt.close(fig2)
#------------------------------------------------------------------
