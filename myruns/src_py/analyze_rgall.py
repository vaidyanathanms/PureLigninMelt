# Analyze and plot Rg/Segmental Rg/Shape factor
# Version: Mar-16-2021
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
# Plot avg rg and rg distribution
print("Analyzing Rg data ...")

for bio_indx in range(len(biom_arr)): # loop in biomass

    biomass = biom_arr[bio_indx]
    res_dir = scr_dir + '/' + inp_type + '/'+ biomass +\
              '/pdi_' + str(pdi_val) + '/results'
    if not os.path.isdir(res_dir):
        os.mkdir(res_dir)
        
    # Write data and then read to concatenate
    rg_fyl = res_dir + '/Rgdata.dat'
    fall_rg = open(rg_fyl,'w')
    fall_rg.write('%s\t%s\t%s\t%s\n'\
                  %('Temperature','<Rg^2>','<Rg^4>','Rg4/Rg2^2'))

    # Define axes labels for plotting Rg distribution
    fig1,ax1 = plt.subplots()
    set_axes(ax1,plt,r'$R_{g}$ (nm)','Probability')
    
    for tval in range(temp_min,temp_max,temp_dt): # loop in temp
        temp_leg  = str(tval)
        rgall = np.zeros(0)
        chain_rg2 = 0; chain_rg4 = 0

        for casenum in range(len(run_arr)): # loop in runarr
            wdir,tdir,fig_dir = ret_temp_dir(scr_dir,inp_type,biomass,\
                                             pdi_val,run_arr[casenum],\
                                             tval,solv_type)
            
            print("Analyzing: ", pdi_val,tval,run_arr[casenum])
            # Check file
            list_of_files = glob.glob(tdir + '/all_results/rg_nptmain_*.xvg')
            if list_of_files == []:
                print("Rg files do not exist for ", tval)
                continue

            if len(list_of_files) != nchains:
                print('ERR: Mismatch in number of analysis file and input',\
                          nchains, len(list_of_files))
                print(list_of_files)
                continue
                
            mon_arr = ret_mons(tdir + '/all_results/chainlist.dat')
            rg_case = res_dir + '/Rg_casenum_' + str(casenum)+ \
                      '_T_' + str(tval) + '.dat'
            fcase_rg = open(rg_case,'w')
            fcase_rg.write('%s\t%s\t%s\t%s\t%s\n' \
                           %('ChainID','Nmons','<Rg^2>','<Rg^4>','Rg4/rg2^2'))

            case_rg2 = 0; case_rg4 = 0
            for fyle in list_of_files: # chain loop
                chid = ret_chid(fyle)
                # Open and parse file
                with open(fyle) as fin:
                    lines = (line.lstrip() for line in fin \
                             if not line.lstrip().startswith('#') and \
                             not line.lstrip().startswith('@'))
                    data  = np.loadtxt(lines)
                    l1 = int(0.1*len(data[:1]))
                    l2 = len(data[:,1])
                    rg2 = np.average(np.square(data[l1:l2,1]))
                    rg4 = np.average(np.power(data[l1:l2,1],4))
                    al =  rg4/(rg2*rg2)
                    fcase_rg.write('%g\t%g\t%g\t%g\t%g\n' %(chid,int(mon_arr[chid]),rg2,rg4,al))

                    case_rg2 += rg2
                    case_rg4 += rg4
                    rgall = np.append(rgall,data[l1:l2,1]) #append all rg-data
                        
            chain_rg2 += case_rg2/nchains
            chain_rg4 += case_rg4/nchains

        chain_rg2 /= len(run_arr)
        chain_rg4 /= len(run_arr)
        alpha = (chain_rg4)/(chain_rg2*chain_rg2)
        fall_rg.write('%g\t%g\t%g\t%g\n' %(tval,chain_rg2,\
                                           chain_rg4,alpha))

        # plot histogram across all cases
        sns.kdeplot(data=np.array(rgall),label=temp_leg,ax=ax1)

        # Save histogram plot
        ax1.legend(loc=0)
        fig1.savefig(fig_dir + '/'+biomass+'_pdi_'+str(pdi_val)+\
                     '_Rgdist.png',dpi=fig1.dpi)
        plt.close(fig1)

    fall_rg.close()
    # plot rg_data
    df=pd.read_table(rg_fyl)
    plot_allrg(df,fig_dir,pdi_val)
#------------------------------------------------------------------
# Plot shape factor and shape factor distribution
print("Analyzing Shape factor data")

for bio_indx in range(len(biom_arr)): # loop in biomass

    biomass = biom_arr[bio_indx]
    res_dir = scr_dir + '/' + inp_type + '/'+ biomass +\
              '/pdi_' + str(pdi_val) + '/results'
    if not os.path.isdir(res_dir):
        os.mkdir(res_dir)
            
    # Write data and then read to concatenate
    sf_fyl = res_dir + '/shapefacdata.dat'
    fall_sf = open(sf_fyl,'w')
    fall_sf.write('%s\t%s\t%s\t%s\t%s\n'\
                  %('Temperature','<lam_1>','<lam_2>','<lam_3>','\kappa'))

    # Define axes labels for plotting Rg distribution
    fig1,ax1 = plt.subplots()
    set_axes(ax1,plt,r'$\kappa$','Probability')

    for tval in range(temp_min,temp_max,temp_dt): # loop in temp
        temp_leg  = str(tval)
        sf_all = np.zeros(0)
        chain_lam1 = 0; chain_lam2 = 0; chain_lam3 = 0; chain_kapa = 0
        
        for casenum in range(len(run_arr)): # loop in runarr
            wdir,tdir,fig_dir = ret_temp_dir(scr_dir,inp_type,biomass,\
                                             pdi_val,run_arr[casenum],\
                                             tval,solv_type)

            print("Analyzing: ", pdi_val,tval,run_arr[casenum])
            # Check file
            list_of_files = glob.glob(tdir + '/all_results/eig_nptmain_*.xvg')
            if list_of_files == []:
                print("Eigenvalue files do not exist for ", tval)
                continue
            
            if len(list_of_files) != nchains:
                print('ERR: Mismatch in number of analysis file and input',\
                      nchains, len(list_of_files))
                print(list_of_files)
                continue
                
            mon_arr = ret_mons(tdir + '/all_results/chainlist.dat')
            sf_case = res_dir + '/sf_casenum_' + str(casenum)+ \
                      '_T_' + str(tval) + '.dat'
            fcase_sf = open(sf_case,'w')
            fcase_sf.write('%s\t%s\t%s\t%s\t%s\t%s\n' \
                           %('ChainID','Nmons','<lam_1>','<lam_2>',\
                             '<lam_3>','\kappa'))
            
            case_lam1 = 0; case_lam2 = 0; case_lam3 = 0; case_kapa = 0
            for fyle in list_of_files:
                chid = ret_chid(fyle)
                # Open and parse file
                with open(fyle) as fin:
                    lines = (line.lstrip() for line in fin \
                             if not line.lstrip().startswith('#') and \
                             not line.lstrip().startswith('@'))
                    data  = np.loadtxt(lines)
                    l1 = int(0.1*len(data[:1]))
                    l2 = len(data[:,1])
                    lam1 = np.average(data[l1:l2,2])
                    lam2 = np.average(data[l1:l2,3])
                    lam2 = np.average(data[l1:l2,4])
                    tr = (lam1+lam2+lam3)/3
                    kappa = 1.5*((lam1-tr)**2 + (lam2-tr)**2 + (lam3-tr)**2)
                    kappa /= (lam1+lam2+lam3)**2
                    fcase_sf.write('%g\t%g\t%g\t%g\t%g\t%g\n' \
                                   %(chid,int(mon_arr[chid]),lam1,lam2,\
                                     lam3,kappa))
                    tr_all = (data[l1:l2,2]+data[l1:l2,3]+data[l1:l2,4])/3
                    sf_all = 1.5*((data[l1:l2,2]-tr_all)**2 + \
                                  (data[l1:l2,3]-tr_all)**2 + \
                                  (data[l2:l2,4]-tr_all)**2))
                    
                    case_lam1 += lam1
                    case_lam2 += lam2
                    case_lam3 += lam3
                    case_kapa += kappa
                    
            chain_lam1 += case_lam1/nchains
            chain_lam2 += case_lam2/nchains
            chain_lam3 += case_lam3/nchains
            chain_kapa += case_kappa/nchains

        chain_lam1 /= len(run_arr)
        chain_lam2 /= len(run_arr)
        chain_lam3 /= len(run_arr)
        chain_kapa /= len(run_arr)
        
        fall_sf.write('%g\t%g\t%g\t%g\t%g\n'\
                      %(tval,chain_lam1,chain_lam2,\
                        chain_lam3,chain_kapa))
        
        # plot histogram across all cases
        sns.kdeplot(data=np.array(sf_all),label=temp_leg,ax=ax1)
        
        # Save histogram plot
        ax1.legend(loc=0)
        fig1.savefig(fig_dir + '/'+biomass+'_pdi_'+str(pdi_val)+\
                     '_sfdist.png',dpi=fig1.dpi)
        plt.close(fig1)
#------------------------------------------------------------------
# Compute Rgscaling
print("Analyzing Segmental Rg data")
for bio_indx in range(len(biom_arr)): # loop in biomass
    biomass = biom_arr[bio_indx]
    res_dir = scr_dir + '/' + inp_type + '/'+ biomass +\
              '/pdi_' + str(pdi_val) + '/results'
    if not os.path.isdir(res_dir):
        os.mkdir(res_dir)
        
    # Write data and then read to concatenate
    nu_fyl = res_dir + '/Rgscaling.dat'
    fall_fit = open(nu_fyl,'w')
    fall_fit.write('%s\t%s\t%s\t%s\n' %('Temp','RunNum','b','nu'))
    
    
    for tval in range(temp_min,temp_max,temp_dt): # loop in temp
        temp_leg  = str(tval)
        
        for casenum in range(len(run_arr)): # loop in runarr
            wdir,tdir,fig_dir = ret_temp_dir(scr_dir,inp_type,biomass,\
                                             pdi_val,run_arr[casenum],\
                                             tval,solv_type)
            
            print("Analyzing: ", pdi_val,tval,run_arr[casenum])
            # Check file
            outfile = tdir + '/all_results/RgvsN.dat'
            if not os.path.exists(outfile):
                print("RgvsN.dat not found for ", tval)
                continue
                    
            df = pd.read_csv(outfile,sep="\s+")
            fits,err = compute_rgscaling(df)
            fall_fit.write('%d\t%d\t%g\t%g\n' \
                           %(tval,run_arr[casenum],fits[0],fits[1]))
    fall_fit.close()

    if len(pdi_arr) > 1: #plot as a function of PDI
        
        fig2, ax2 = plt.subplots()
        set_axes(ax2,plt,r'Temperature ($K$)',r'$\nu$') 

        fig3, ax3 = plt.subplots()
        set_axes(ax2,plt,r'Temperature ($K$)',r'$C_{l} (nm)$') 

        for pdiindx in range(len(pdi_arr)):
            pdiplt = pdi_arr[pdiindx]
            pdileg = 'PDI: ' + str(pdiplt)
            res_dir = scr_dir + '/' + inp_type + '/'+ biomass +\
                      '/pdi_' + str(pdiplt) + '/results'

            if not os.path.isdir(res_dir):
                print(res_dir, "does not exist")
                continue

            else:
                rg_fyl = res_dir + '/Rgscaling.dat'
                df=pd.read_csv(rg_fyl,"\s+")
                ax2.scatter(x=df['T'],y=df['nu'],marker=mrk_arr[pdiindx],\
                            color=clr_arr[pdiindx])
                ax3.scatter(x=df['T'],y=df['b'],marker=mrk_arr[pdiindx],\
                            color=clr_arr[pdiindx])

        fig2.savefig(fig_dir + '/'+'nu_all.png',dpi=fig2.dpi)
        plt.close(fig2)
        fig3.savefig(fig_dir + '/'+'cl_all.png',dpi=fig3.dpi)
        plt.close(fig3)
#------------------------------------------------------------------
