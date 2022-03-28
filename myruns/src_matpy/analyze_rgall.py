# Analyze and plot Rg
# Version: Mar-27-2022
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
temp_max  = 301 # Maximum temperature
temp_dt   = 10  # Temperature dt
pdi_arr   = [1.0,1.8,3.0,'expts']
mark_arr  = ['o','d','s']
nchains   = 20
anadir_head = 'all_radgyr' # change this for different analysis
start_frac  = 0.3 # starting point for averaging
#------------------------------------------------------------------

# Global arrays
denarr = np.arange(temp_min,temp_max,5*temp_dt) #plotting time data
#------------------------------------------------------------------

# Generate output directories
anaout_dir = anaout_dir + '/rg_results' # result_outputs
if not os.path.isdir(anaout_dir):
    os.mkdir(anaout_dir)
#------------------------------------------------------------------

#----------------Begin Analyzing Rg2/Rg4 data -----------------------------
print("Analyzing Rg data ...")
for pdi_val in pdi_arr:

    # Set/Check directories
    simout_dir = all_dir + '/results_pdi_' + str(pdi_val)
    if not os.path.isdir(simout_dir):
        print("ERR: " +  simout_dir + " does not exist")       
        continue
    
    # Output for each case
    fall_out = open(anaout_dir +'/allRg2data_'+str(pdi_val)+'.dat','w')
    fall_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %('Temp','Rg2_R1','Rg2_R2',\
                                                'Rg2_R3','Rg2_R4','TotCase'))

    fall2_out = open(anaout_dir +'/allRg4data_'+str(pdi_val)+'.dat','w')
    fall2_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %('Temp','Rg4_R1','Rg4_R2',\
                                                 'Rg4_R3','Rg4_R4','TotCase'))


    fall3_out = open(anaout_dir +'/allalpha_'+str(pdi_val)+'.dat','w')
    fall3_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %('Temp','Al_R1','Al_R2',\
                                                 'Al_R3','Al_R4','TotCase'))

    # Averaged outputs
    fc_rg = open(anaout_dir +'/Rgdata_'+str(pdi_val)+'.dat','w')
    fc_rg.write('%s\t%s\t%s\t%s\t%s\n' \
                %('Temp','<Rg^2>','<Rg^4>','Rg4/Rg2^2','Ncases'))

    # Define axes labels for distribution
    fig1,ax1 = plt.subplots()
    set_axes(ax1,plt,r'$R_{g}$ (nm)','Probability')

    pdiflag = -1                # to account for empty directories
    
    # Temperature loops and averaging
    for tval in range(temp_min,temp_max,temp_dt): # loop in temp
        temp_leg  = str(tval)
        rgall = []; pdiflag = 1
        chain_rg2 = 0; chain_rg4 = 0; ncases_pertemp = 0; alpha = 0
        fall_out.write('%g\t' %(tval)); fall2_out.write('%g\t' %(tval))
        fall3_out.write('%g\t' %(tval))
        
        for casenum in range(len(run_arr)): # loop in runarr
            wdir = simout_dir + '/run_' + str(run_arr[casenum]) +\
                '/T_' + str(tval) + '/' + anadir_head
            if not os.path.isdir(wdir):
                fall_out.write('%s\t' %('N/A')); fall2_out.write('%s\t' %('N/A'))
                fall3_out.write('%s\t' %('N/A'))
                print("ERR: " + wdir + " does not exist")
                continue

            
            print("Analyzing: ", pdi_val,tval,run_arr[casenum])
            
            # Check file(s)
            list_of_files = glob.glob(wdir + '/rg_nptmain_*.xvg')
            if list_of_files == []:
                fall_out.write('%s\t' %('N/A')); fall2_out.write('%s\t' %('N/A'))
                fall3_out.write('%s\t' %('N/A'))
                print("Rg files do not exist for ", tval)
                continue

            if len(list_of_files) != nchains:
                fall_out.write('%s\t' %('DiffNchains'))
                fall2_out.write('%s\t' %('DiffNchains'))
                fall3_out.write('%s\t' %('DiffNchains'))
                print('ERR: Mismatch in number of analysis file and input',\
                          nchains, len(list_of_files))
                print(list_of_files)
                continue

            if not os.path.exists(wdir + '/chainlist.dat'):
                fall_out.write('%s\t' %('Nochainlist'))
                fall2_out.write('%s\t' %('Nochainlist'))
                fall3_out.write('%s\t' %('Nochainlist'))
                print('ERR: chainlist.dat not found')
                continue
            
            mon_arr = ret_mons(wdir + '/chainlist.dat')
            fcase_rg = open(anaout_dir + '/Rg_casenum_' + str(run_arr[casenum])+ \
                            '_T_' + str(tval) + '.dat','w')
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
                    data  = np.loadtxt(lines) #end with open(fyle)
                    
                l1 = int(start_frac*len(data[:1]))
                l2 = len(data[:,1])
                rg2 = np.average(np.square(data[l1:l2,1]))
                rg4 = np.average(np.power(data[l1:l2,1],4))
                fcase_rg.write('%g\t%g\t%g\t%g\t%g\n'\
                               %(chid,int(mon_arr[chid]),rg2,rg4,rg4/rg2))
                 
                case_rg2 += rg2
                case_rg4 += rg4

                #append new rg-data
                rgall = np.append(rgall,data[l1:l2,1]) #end for fyle in list_of_files

            fcase_rg.close() # close writing for each case 
            chain_rg2 += case_rg2/nchains
            chain_rg4 += case_rg4/nchains
            al         = (nchains*case_rg4)/(case_rg2*case_rg2)
            alpha += al 
            ncases_pertemp += 1 

            # Write average rg2/rg4/alpha for each case
            fall_out.write('%g\t' %(case_rg2/nchains))
            fall2_out.write('%g\t' %(case_rg4/nchains)) 
            fall3_out.write('%g\t' %(al))
            # end for casenum in len(range(run_arr))
        
        # Do NOT continue if zero cases are found
        if ncases_pertemp == 0:
            fall_out.write('%g\n' %(ncases_pertemp))
            fall2_out.write('%g\t' %(ncases_pertemp))
            fall3_out.write('%g\t' %(ncases_pertemp))             
            continue

        # Write averages over different cases for each temp
        chain_rg2 /= ncases_pertemp
        chain_rg4 /= ncases_pertemp
        alpha     /= ncases_pertemp
        fc_rg.write('%g\t%g\t%g\t%g\t%d\n' %(tval,chain_rg2,\
                                             chain_rg4,alpha,ncases_pertemp))
        fall_out.write('%g\n' %(ncases_pertemp))
        fall2_out.write('%g\n' %(ncases_pertemp))
        fall3_out.write('%g\n' %(ncases_pertemp)) 

        # plot Rg-distribution if tval in denarr
        if tval in denarr: # end of tval loop
            sns.kdeplot(data=np.array(rgall),label=temp_leg,ax=ax1)
            ax1.legend(loc=0) 
   
    # Close temp files
    fc_rg.close();fall_out.close(); fall2_out.close(); fall3_out.close()

    # Do NOT continue if no PDI-temps are found
    if pdiflag == -1: 
        continue

    # Save distribution for select temperatures
    fig1.savefig(figout_dir+'/'+'Rgdist_'+str(pdi_val)+'.png',dpi=fig1.dpi)
    fig1.savefig(figout_dir+'/'+'Rgdist_'+str(pdi_val)+'.eps',format='eps')
    plt.close(fig1)

    # Plot Rg-temperature data for each PDI
    df=pd.read_table(anaout_dir +'/Rgdata_'+str(pdi_val)+'.dat')
    plot_allrg(df,figout_dir,pdi_val) #end of PDI loop
    
#------- Plot Rg2-Temp data for all PDI values together--------------------
print("Plotting all Rg data as a function of PDI..")
fig2, ax2 = plt.subplots()
set_axes(ax2,plt,r'Temperature ($K$)',r'\langle R_{g}^{2} \rangle ($nm^{2}$)')

for pdi_val in pdi_arr:
    if pdi_val == 'expts':
        pdileg = 'PDI: Experimental Distribution'
    else:
        pdileg = 'PDI: ' + str(pdi_val)
    fname = '/Rgdata_'+str(pdi_val)+'.dat'
    if not os.path.exists(anaout_dir + '/' + fname):
        print('ERR: '+fname+' does not exist in ' + anaout_dir)
        continue
    
    df=pd.read_table(anaout_dir + '/' + fname)
    print('Plotting', pdi_val)
    ax2.scatter(x=df['Temp'],y=df['<Rg^2>'],label=pdileg)
    
fig2.savefig(figout_dir + '/'+'rg2_allpdi.png',dpi=fig2.dpi)
fig2.savefig(figout_dir + '/'+'rg2_allpdi.eps',format='eps')
plt.close(fig2)
#----------------------End of Rg2 Analysis------------------------------
