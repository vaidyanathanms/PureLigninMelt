# Analyze and plot Shape factor
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

#----- Begin shape factor analysis --------------------------------
print("Analyzing Shape factor data")

for pdi_val in pdi_arr:

    # Set/Check directories
    simout_dir = all_dir + '/results_pdi_' + str(pdi_val)
    if not os.path.isdir(simout_dir):
        print("ERR: " +  simout_dir + " does not exist")       
        continue

    # Output for each case
    fall_out = open(anaout_dir +'/allRg2data_'+str(pdi_val)+'.dat','w')
    fall_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %('Temp','kap_R1','kap_R2',\
                                                'kap_R3','kap_R4','TotCase'))


    # Averaged outputs
    sf_fyl  = anaout_dir + '/shapefacdata_'+str(pdi_val)+'.dat' 
    fall_sf = open(sf_fyl,'w')
    fall_sf.write('%s\t%s\t%s\t%s\t%s\t%s\n'\
                  %('Temp','<lam_1>','<lam_2>','<lam_3>','<\kappa>',\
                    'Ncases'))

    # Define axes labels for shape factor distribution
    fig1,ax1 = plt.subplots()
    set_axes(ax1,plt,r'$\kappa$','Probability')

    # Temperature loops and averaging
    for tval in range(temp_min,temp_max,temp_dt): # loop in temp
        temp_leg  = str(tval)
        sf_all = []
        avgd_lam1 = 0; avgd_lam2 = 0; avgd_lam3 = 0; avgd_kappa = 0
        ncases_pertemp = 0
        fall_out.write('%g\t' %(tval))
        
        for casenum in range(len(run_arr)): # loop in runarr
            wdir = simout_dir + '/run_' + str(run_arr[casenum]) +\
                '/T_' + str(tval) + '/' + anadir_head
            if not os.path.isdir(wdir):
                fall_out.write('%s\t' %('N/A'))
                print("ERR: " + wdir + " does not exist")
                continue


            print("Analyzing: ", pdi_val,tval,run_arr[casenum])

            # Check file(s)
            list_of_files = glob.glob(wdir + '/eig_nptmain_*.xvg')
            if list_of_files == []:
                fall_out.write('%s\t' %('N/A'));
                print("Eigenvalue files do not exist for ", tval)
                continue
            
            if len(list_of_files) != nchains:
                fall_out.write('%s\t' %('DiffNchains'));
                print('ERR: Mismatch in number of analysis file and input',\
                      nchains, len(list_of_files))
                print(list_of_files)
                continue

            if not os.path.exists(wdir + '/chainlist.dat'):
                fall_out.write('%s\t' %('Nochainlist'))
                print('ERR: chainlist.dat not found')
                continue

            mon_arr = ret_mons(wdir + '/chainlist.dat')
            sf_case = anaout_dir + '/shape_casenum_' + str(run_arr[casenum])+ \
                      '_T_' + str(tval) + '.dat'
            fcase_sf = open(sf_case,'w')
            fcase_sf.write('%s\t%s\t%s\t%s\t%s\t%s\n' \
                           %('ChainID','Nmons','<lam_1>','<lam_2>',\
                             '<lam_3>','\kappa'))
            
            case_lam1 = 0; case_lam2 = 0; case_lam3 = 0; case_kapa = 0

            for fyle in list_of_files: # chain loop
                chid = ret_chid(fyle)
                # Open and parse file
                with open(fyle) as fin:
                    lines = (line.lstrip() for line in fin \
                             if not line.lstrip().startswith('#') and \
                             not line.lstrip().startswith('@'))
                    data  = np.loadtxt(lines) #end with open(fyle)

                l1 = int(0.1*len(data[:1]))
                l2 = len(data[:,1])
                lam1 = np.average(data[l1:l2,2])
                lam2 = np.average(data[l1:l2,3])
                lam3 = np.average(data[l1:l2,4])
                tr   = (lam1+lam2+lam3)/3
                kappa = 1.5*((lam1-tr)**2 + (lam2-tr)**2 + (lam3-tr)**2)
                kappa /= (lam1+lam2+lam3)**2
                fcase_sf.write('%g\t%g\t%g\t%g\t%g\t%g\n' \
                               %(chid,int(mon_arr[chid]),lam1,lam2,\
                                 lam3,kappa))

                case_lam1 += lam1
                case_lam2 += lam2
                case_lam3 += lam3
                case_kapa += kappa 

                # append shape factor data
                sf_all = np.append(sf_all,kappa)
                                   
                #end for fyle in list_of_files                 
                
            fcase_sf.close() #close writing for each case    
            avgd_lam1 += case_lam1/nchains
            avgd_lam2 += case_lam2/nchains
            avgd_lam3 += case_lam3/nchains
            avgd_kappa += case_kapa/nchains
            ncases_pertemp += 1
            
            fall_out.write('%g\t' %(case_kapa/nchains)) #write each case kappa
            # end for casenum in len(range(run_arr))

        # Do NOT continue if zero cases are found
        if ncases_pertemp == 0:
            fall_out.write('%g\n' %(ncases_pertemp))
            continue

        # Write averages over different cases for each temp
        avgd_lam1 /= ncases_pertemp
        avgd_lam2 /= ncases_pertemp
        avgd_lam3 /= ncases_pertemp
        avgd_kappa /= ncases_pertemp
        fall_out.write('%g\n' %(ncases_pertemp))
        fall_sf.write('%g\t%g\t%g\t%g\t%g\t%d\n' %(tval,avgd_lam1,avgd_lam2,\
                                                   avgd_lam3,avgd_kappa,\
                                                   ncases_pertemp))

        # plot kappa-distribution across all cases
        if tval in denarr: # end of tval loop
            sns.kdeplot(data=np.array(sf_all),label=temp_leg,ax=ax1)
            ax1.legend(loc=0)

    # Close temp files
    fall_out.close(); fall_sf.close()
    
    # Do NOT continue if no PDI-temps are found
    if len(sf_all) == 0: #to account for boolean values
        continue

    # Save distribution for select temperatures
    ax1.legend(loc=0)
    fig1.savefig(figout_dir+'/'+'kapdist_'+str(pdi_val)+'.png',dpi=fig1.dpi)
    fig1.savefig(figout_dir+'/'+'kapdist_'+str(pdi_val)+'.eps',format='eps')
    plt.close(fig1)

    # Plot Rg-temperature data for each PDI
    df=pd.read_table(sf_fyl)
    plot_allshape(df,figout_dir,pdi_val) #end of PDI loop

#------- Plot Kappa-Temp data for all PDI values together--------------------
print("Plotting all shape data as a function of PDI..")
fig2, ax2 = plt.subplots()
set_axes(ax2,plt,r'Temperature ($K$)',r'\langle \kappa \rangle')

for pdi_val in pdi_arr:
    if pdi_val == 'expts':
        pdileg = 'PDI: Experimental Distribution'
    else:
        pdileg = 'PDI: ' + str(pdi_val)
    fname = '/shapefacdata_'+str(pdi_val)+'.dat' 
    if not os.path.exists(anaout_dir + '/' + fname):
        print('ERR: '+fname+' does not exist in ' + anaout_dir)
        continue
    
    df=pd.read_table(anaout_dir + '/' + fname)
    print('Plotting', pdi_val)
    ax2.scatter(x=df['Temp'],y=df['<\kappa>'],label=pdileg)
    
fig2.savefig(figout_dir + '/'+'kap_allpdi.png',dpi=fig2.dpi)
fig2.savefig(figout_dir + '/'+'kap_allpdi.eps',format='eps')
plt.close(fig2)
#----------------------End of kappa Analysis------------------------------

