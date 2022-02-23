# Generic plot
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

# Directory paths
main_dir = os.getcwd() # current dir
scr_dir  = '/lustre/or-scratch/cades-bsd/v0e' # scratch dir
scr_dir  = scr_dir + '/Glassy_lignin'
sh_dir   = '../src_gmx/sh_files'
if not os.path.isdir(scr_dir):
    print("FATAL ERROR: ", scr_dir, " not found")
    exit("Check scratch directory path")
#------------------------------------------------------------------

# Plot keys 0,1,2 (0-None,1-separate,2-together)
rg_cal = 0 # avg rg and rg distribution
tg_cal = 0 # plot specific volume for tg
rdf_pol = 0 # plot RDFs
rg_scaling  = 0 # compute and plot rg scaling
den_cal = 1 # plot densities
sf_cal = 0 # shape factor
#------------------------------------------------------------------

# Color/line data; figure defaults
orange = '#FFA500'; dark_g = '#006400'; brown = '#8B4513'
clr_arr = [dark_g,orange,brown,'b','k','m']
mrk_arr = ['o','d','s']
lne_arr = ['-','--']
plt.rc('legend',fontsize=16) # fontsize of the legends
plt.rcParams.update({'font.size': 16}) # other fonts
plt.rcParams.update({'figure.autolayout': True})
#plt.rcParams['font.family'] = 'serif'
#plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
#plt.rcParams.update({'font.family': 'Times New Roman'})
#------------------------------------------------------------------

# Input data
inp_type  = 'melts' # melts, solvents, cosolvents
disperse  = 'poly' # mono/poly; only for melts
biom_arr  = ['WT'] # biomass type arr
otyp_arr  = ['None']  # solvent arr for solvents/cosolvents
otyp_leg  = ['None']  # legend for solvents/cosolvents
solv_type = 'None'
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

# Compute densities for all systems
if den_cal != 0:
    print("Computing densities")
    if not os.path.isdir(sh_dir):
        raise RuntimeError(sh_dir + " does not exist!")
    for bio_indx in range(len(biom_arr)): # loop in biomass
        biomass = biom_arr[bio_indx]
        for casenum in range(len(run_arr)): # loop in runarr
            wdir,tdir,fig_dir = ret_temp_dir(scr_dir,inp_type,biomass,\
                                             pdi_val,run_arr[casenum],\
                                             600,solv_type)
            print(wdir)
            gencpy(sh_dir,wdir,'comp_dens_pyinp.sh')
            submit_dens(os.getcwd(),wdir,'comp_dens_pyinp.sh',\
                        temp_min,temp_max,temp_dt)
#------------------------------------------------------------------

# Plot glass transition for all systems
if tg_cal != 0:
    print("Analyzing Tg data")

    for bio_indx in range(len(biom_arr)): # loop in biomass

        biomass = biom_arr[bio_indx]
        if tg_cal == 2: # plot all averages together
            res_dir = scr_dir + '/' + inp_type + '/'+ biomass +\
                      '/pdi_' + str(pdi_val) + '/results'
            if not os.path.isdir(res_dir):
                os.mkdir(res_dir)

            # Write data and then read to concatenate
            tg_fyl = res_dir + '/tgdata.dat'
            fc_tg = open(tg_fyl,'w')
            fc_tg.write('%s\t%s\t%s\n' \
                        %('Biomass','Temp','SV_NPT'))


        tg_avg = np.zeros(0)

        # Define axes labels for distribution
        fig1,ax1 = plt.subplots()
        set_axes(ax1,plt,r'Time (ps)',r'Density ($kg$/$m^3$)')

        ytemp_avg = np.zeros(0)

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
                halflen = int(0.5*len(yall))
                yavg += np.average(data[halflen:len(yall)-1,1])

                if casenum == 0 and tval in denarr:
                    ax1.plot(data[:,0],data[:,1],label=temp_leg)

            # append avg to yavg
            ytemp_avg = np.append(ytemp_avg,yavg/len(run_arr))

        # Save density plots
        ax1.legend(loc=0)
        fig1.savefig(fig_dir + '/'+'dens_'+str(pdi_val)+'.png',\
                     dpi=fig1.dpi)
        plt.close(fig1)

        if tg_cal == 2:
            tval = 300
            for indx in range(len(ytemp_avg)): 
                temp_tot  = ytemp_avg[indx]
                if temp_tot == 0:
                    continue
                sval = 1000.0/temp_tot
                fc_tg.write('%s\t%g\t%g\n' %(biomass,tval,\
                                             sval))
                tval += temp_dt

    if tg_cal == 2:
        fc_tg.close()
        df=pd.read_table(tg_fyl)
        figa, axa = plot_tg(df,fig_dir,pdi_val) 
        tgval = compute_tg(df,axa) #compute tgval
        figa.savefig(fig_dir + '/'+'sv2_' + str(pdi_val) + '.png',\
                     dpi=figa.dpi)
        plt.close(figa)

    if len(pdi_arr) > 1: #plot as a function of PDI

        fig2, ax2 = plt.subplots()
        set_axes(ax2,plt,r'Temperature ($K$)',r'Specific Volume ($cm^3/g$)')

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
# Plot avg rg and rg distribution
if rg_cal != 0:
    print("Analyzing Rg data")

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
                for fyle in list_of_files:
                    chid = ret_chid(fyle)
                    # Open and parse file
                    with open(fyle) as fin:
                        lines = (line.lstrip() for line in fin \
                                 if not line.lstrip().startswith('#') and \
                                 not line.lstrip().startswith('@'))
                        data  = np.loadtxt(lines)
                        l1 = int(0.2*len(data[:1]))
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
            if rg_cal == 2:
                fall_rg.write('%g\t%g\t%g\t%g\n' %(tval,chain_rg2,\
                                                   chain_rg4,alpha))

            # plot histogram across all cases
            sns.kdeplot(data=np.array(rgall),label=temp_leg,ax=ax1)

        # Save histogram plot
        ax1.legend(loc=0)
        fig1.savefig(fig_dir + '/'+biomass+'_pdi_'+str(pdi_val)+\
                     '_Rgdist.png',dpi=fig1.dpi)
        plt.close(fig1)

    if rg_cal == 2:
        fall_rg.close()
        df=pd.read_table(rg_fyl)
        plot_allrg(df,fig_dir,pdi_val)
#------------------------------------------------------------------
# Plot shape factor and shape factor distribution
if sf_cal != 0:
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
                        l1 = int(0.2*len(data[:1]))
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

            if sf_cal == 2:
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
if rg_scaling != 0:
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
        fall_fit.write('%s\t%s\t%s\n' \
                   %('T','b','nu'))


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
                fall_fit.write('%d\t%g\t%g\n' \
                               %(tval,fits[0],fits[1]))
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
