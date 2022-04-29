# Analyze and plot MSDs
# Version: Apr-26-2022
#------------------------------------------------------------------

# Import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
import math
from plt_aux import *
import compute_props as cprop
from plt_inps import *
#------------------------------------------------------------------

# Input data for analysis
run_arr   = [1,2,3,4] # run numbers for a given biomass
temp_min  = 250 # Minimum temperature
temp_max  = 501 # Maximum temperature
temp_dt   = 50  # Temperature dt
pdi_arr   = [1.0,1.8,3.0,'3.7','expts']
nchains   = 20
anadir_head = 'all_msd' # change this for different analysis
tref = 100000 # reference time in ps 
#------------------------------------------------------------------

# Generate output directories
anaout_dir = anaout_dir + '/msd_results' # result_outputs
if not os.path.isdir(anaout_dir):
    os.mkdir(anaout_dir)
#------------------------------------------------------------------

#----------------Begin Analyzing MSD data--------------------------
print("Analyzing MSD data ...")
for pdi_val in pdi_arr:

    # Set/Check directories
    simout_dir = all_dir + '/results_pdi_' + str(pdi_val)
    if not os.path.isdir(simout_dir):
        print("ERR: " +  simout_dir + " does not exist")       
        continue
    
    # Temperature loops and averaging
    for tval in range(temp_min,temp_max,temp_dt): # loop in temp
        temp_leg  = str(tval)
        pdiflag = 1;  msdtime = 0

        # Output at a given temperature for all cases
        fall_out = open(anaout_dir +'/allmsddata_T_'+str(tval)+'_pdi_'\
                        + str(pdi_val)+'.dat','w')
        fall_out.write('%s\t%g\t%s\n' %('MSD @t = ', tref, ' ps'))
        fall_out.write('%s\t%s\t%s\t%s\n' %('CaseNum','ChID','NMons','MSD'))

        for casenum in range(len(run_arr)): # loop in runarr
            wdir = simout_dir + '/run_' + str(run_arr[casenum]) +\
                '/T_' + str(tval) + '/' + anadir_head
            
            print("Analyzing: ", pdi_val,tval,run_arr[casenum])
            
            fcase_msd = open(anaout_dir + '/msd_pdi_' + str(pdi_val) + '_casenum_'\
                            + str(run_arr[casenum])+ '_T_' + str(tval)
                            + '.dat','w')
            fcase_msd.write('%s\t%g\t%s\n' %('MSD @t = ', tref, ' ps'))
            fcase_msd.write('%s\t%s\t%s\n' %('ChainID','Nmons','MSD'))

            # Check file(s)
            if not os.path.isdir(wdir):
                fall_out.write('%g\t%s\t%s\t%s\n' %(casenum,'N/A','N/A','N/A'))
                fcase_msd.write('%s\n' %('No directory found'))
                fcase_msd.close()
                print("ERR: " + wdir + " does not exist")
                continue

            msd_list_of_files = glob.glob(wdir + '/msd_nptmain_*.xvg')
            if msd_list_of_files == []:
                fall_out.write('%g\t%s\t%s\t%s\n' %(run_arr[casenum],\
                                                    'N/A','N/A','N/A'))
                fcase_msd.write('%s\n' %('No chains found'))
                fcase_msd.close()
                print("MSD files do not exist for ", tval)
                continue

            if len(msd_list_of_files) != nchains:
                fall_out.write('%g\t%s\t%s\t%s\n' %(run_arr[casenum],\
                                                    'N/A','N/A','N/A'))
                fcase_msd.write('%s\n' %('DiffNchains'))
                fcase_msd.close()
                print('ERR: Mismatch in number of analysis file and input',\
                      nchains, len(msd_list_of_files))
                continue

            if not os.path.exists(wdir + '/chainlist.dat'):
                fall_out.write('%g\t%s\t%s\t%s\n' %(run_arr[casenum],\
                                                    'N/A','N/A','N/A'))
                fcase_msd.write('%s\n' %('Nochainlist'))
                fcase_msd.close()
                print('ERR: chainlist.dat not found')
                continue

            mon_arr = ret_mons(wdir + '/chainlist.dat')

            for fyle in msd_list_of_files: # msd chain loop
                chid = ret_chid(fyle) 
                # Open and parse msd file
                with open(fyle) as fin:
                    lines = (line.lstrip() for line in fin \
                             if not line.lstrip().startswith('#') and \
                             not line.lstrip().startswith('@'))
                    data  = np.loadtxt(lines) #end with open(fyle)
                    
                msd = cprop.msd_tref(data,tref)                   
                fcase_msd.write('%g\t%g\t%g\n' %(chid,mon_arr[chid],msd))
                fall_out.write('%g\t%g\t%g\t%g\n' %(run_arr[casenum],chid,\
                                                    mon_arr[chid],msd))
                 
            # end for fyle in msd_list_of_files
            fcase_msd.close() # close writing for each case
            
        #end casenum in range(len(run_arr))
        fall_out.close()
        
    #end tval in range(temp_min,temp_max,temp_dt)
#end pdival loop
#------- End MSD Analysi-------------------------------------------------
