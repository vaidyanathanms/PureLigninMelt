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
temp_max  = 501 # Maximum temperature
temp_dt   = 50  # Temperature dt
pdi_arr   = [1.0,1.8,3.0,3.7,'expts']
mark_arr  = ['o','d','s']
nchains   = 20
anadir_head = 'all_radgyr' # change this for different analysis
start_frac  = 0.3 # starting point for averaging
#------------------------------------------------------------------

# Global arrays
denarr = np.arange(temp_min,temp_max,5*temp_dt) #plotting time data
#------------------------------------------------------------------

# Generate output directories
anaout_dir = anaout_dir + '/shape_results' # result_outputs
if not os.path.isdir(anaout_dir):
    os.mkdir(anaout_dir)
#------------------------------------------------------------------

#-----Begin shape factor analysis----------------------------------
print("Analyzing Shape factor data")

for pdi_val in pdi_arr:

    # Set/Check directories
    simout_dir = all_dir + '/results_pdi_' + str(pdi_val)
    if not os.path.isdir(simout_dir):
        print("ERR: " +  simout_dir + " does not exist")       
        continue

    pdiflag = -1                # to account for empty directories
    
    # Define axes labels for shape factor distribution
    fig1,ax1 = plt.subplots()
    set_axes(ax1,plt,r'$\kappa$',r'Probability')

    # Temperature loops and averaging
    for tval in range(temp_min,temp_max,temp_dt): # loop in temp
        temp_leg  = str(tval)
        sf_all = []; pdiflag = 1

        # Output at a given temperature for all cases
        fall_out = open(anaout_dir + '/shape_T_' + str(tval)+ '_pdi_'\
                        + str(pdi_val) + '.dat','w')
        fall_out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                       %('CaseNum','ChainID','Nmons','<lam_1>',\
                         '<lam_2>','<lam_3>','\kappa'))
        
        for casenum in range(len(run_arr)): # loop in runarr
            wdir = simout_dir + '/run_' + str(run_arr[casenum]) +\
                '/T_' + str(tval) + '/' + anadir_head
            if not os.path.isdir(wdir):
                fall_out.write('%s\t' %('N/A'))
                print("ERR: " + wdir + " does not exist")
                continue

            print("Analyzing: ", pdi_val,tval,run_arr[casenum])

            fcase_sf = open(anaout_dir + '/shapefac_pdi_' + \
                            str(pdi_val) + '_casenum_' + \
                            str(run_arr[casenum]) + '_T_' + str(tval) \
                            + '.dat','w')
            fcase_sf.write('%s\t%s\t%s\t%s\t%s\t%s\n' \
                           %('ChainID','Nmons','<lam_1>','<lam_2>',\
                             '<lam_3>','\kappa'))

            # Check file(s)
            list_of_files = glob.glob(wdir + '/eig_nptmain_*.xvg')
            if list_of_files == []:
                fall_out.write('%g\t%s\t%s\t%s\t%s\t%s\t%s\n'
                               %(run_arr[casenum],'N/A','N/A','N/A',\
                                 'N/A','N/A','N/A'))
                fcase_sf.write('%s\n' %('No eigenvalue files found'))
                fcase_sf.close()
                print("Eigenvalue files do not exist for ", tval)
                continue
            
            if len(list_of_files) != nchains:
                fall_out.write('%g\t%s\t%s\t%s\t%s\t%s\t%s\n'
                               %(run_arr[casenum],'N/A','N/A','N/A',\
                                 'N/A','N/A','N/A'))
                fcase_sf.write('%s\n' %('Mismatch in number of chains'))
                fcase_sf.close()
                print('ERR: Mismatch in number of chains in file and inp',\
                      nchains, len(list_of_files))
                print(list_of_files)
                continue

            if not os.path.exists(wdir + '/chainlist.dat'):
                fall_out.write('%g\t%s\t%s\t%s\t%s\t%s\t%s\n'
                               %(run_arr[casenum],'N/A','N/A','N/A',\
                                 'N/A','N/A','N/A'))
                fcase_sf.write('%s\n' %('chainlist.dat not found'))
                fcase_sf.close()
                print('ERR: chainlist.dat not found')
                continue

            mon_arr = ret_mons(wdir + '/chainlist.dat')           
            case_lam1 = 0; case_lam2 = 0; case_lam3 = 0; case_kapa = 0

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
                lam1 = np.average(data[l1:l2,2])
                lam2 = np.average(data[l1:l2,3])
                lam3 = np.average(data[l1:l2,4])
                tr   = (lam1+lam2+lam3)/3
                kappa = 1.5*((lam1-tr)**2 + (lam2-tr)**2 + (lam3-tr)**2)
                kappa /= (lam1+lam2+lam3)**2
                fcase_sf.write('%g\t%g\t%g\t%g\t%g\t%g\n' \
                               %(chid,int(mon_arr[chid]),\
                                 lam1,lam2,lam3,kappa))
                fall_out.write('%g\t%g\t%g\t%g\t%g\t%g\t%g\n' \
                               %(casenum,chid,int(mon_arr[chid]),\
                                 lam1,lam2,lam3,kappa))

            #end for fyle in list_of_files                 
            fcase_sf.close()
            
        # end for casenum in len(range(run_arr))
        fall_out.close()
        
    #end tval in range(temp_min,temp_max,temp_dt)
#end pdival loop
#----------------------End of kappa Analysis------------------------------

