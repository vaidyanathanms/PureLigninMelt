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
temp_max  = 501 # Maximum temperature
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

#----------------Begin Analyzing Segmental Rg------------------------

print("Analyzing Segmental Rg data...")
for pdi_val in pdi_arr:

    # Set/Check directories
    simout_dir = all_dir + '/results_pdi_' + str(pdi_val)
    if not os.path.isdir(simout_dir):
        print("ERR: " +  simout_dir + " does not exist")       
        continue
        
    # Output for each case
    fall_out = open(anaout_dir +'/allsegbdata_'+str(pdi_val)+'.dat','w')
    fall_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %('Temp','b_R1','b_R2',\
                                                'b_R3','b_R4','TotCase'))

    fall2_out = open(anaout_dir +'/allnudata_'+str(pdi_val)+'.dat','w')
    fall2_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %('Temp','nu_R1','nu_R2',\
                                                'nu_R3','nu_R4','TotCase'))


    # Averaged outputs
    fall_fit = open(anaout_dir +'/Rgscaling_'+str(pdi_val)+'.dat','w')
    fall_fit.write('%s\t%s\t%s\t%s\n' %('Temp','b','nu','NCases'))

    # Temperature loops and averaging
    for tval in range(temp_min,temp_max,temp_dt): # loop in temp
        temp_leg  = str(tval)
        nu_avg = 0; b_avg = 0; ncases_pertemp = 0
        fall_out.write('%g\t' %(tval)); fall2_out.write('%g\t' %(tval))
        
        for casenum in range(len(run_arr)): # loop in runarr
            wdir = simout_dir + '/run_' + str(run_arr[casenum]) +\
                '/T_' + str(tval) + '/' + anadir_head
            if not os.path.isdir(wdir):
                fall_out.write('%s\t' %('N/A')); fall2_out.write('%s\t' %('N/A'))
                print("ERR: " + wdir + " does not exist")
                continue
            
            print("Analyzing: ", pdi_val,tval,run_arr[casenum])

            # Check file(s)
            list_of_files = wdir + '/RgvsN.dat'
            if not os.path.exists(outfile):
                fall_out.write('%s\t' %('N/A')); fall2_out.write('%s\t' %('N/A'))
                print("RgvsN.dat not found for ", tval)
                continue
                    
            df = pd.read_csv(outfile,sep="\s+")
            fits,err = compute_rgscaling(df)
            fall_out.write('%g\t' %(fits[0])); fall2_out.write('%g\t' %(fits[1]))
            b_avg += fits[0]; nu_avg += fits[1] # end of casenum loop
            ncases_pertemp += 1
            
        # Do NOT continue if zero cases are found
        if ncases_pertemp == 0:
            fall_out.write('%g\n' %(ncases_pertemp))
            fall2_out.write('%g\n' %(ncases_pertemp))
            continue

        # Average and write for each temp
        b_avg /= ncases_pertemp
        nu_avg /= ncases_pertemp
        fall_out.write('%g\n' %(ncases_pertemp))
        fall2_out.write('%g\n' %(ncases_pertemp))
        fall_fit.write('%g\t%g\t%g\t%g\n' %(tval,b_avg,nu_avg,ncases_pertemp))
        # end of tval loop

    # Close temp files
    fall_fit.close(); fall_out.close(); fall_out2.close()
    
    # Do NOT continue if no PDI-temps are found
    if len(rgall) == 0: #to account for boolean values
        continue
        

    # Plot b-Temp data and nu-temp data
    df=pd.read_table(anaout_dir +'/Rgscaling_'+str(pdi_val)+'.dat')
    plot_rgscaling(df,figout_dir,pdi_val) #end of PDI loop

#------- Plot b/nu-Temp data for all PDI values together--------------------

fig2, ax2 = plt.subplots()
set_axes(ax2,plt,r'Temperature ($K$)',r'$\nu$') 

fig3, ax3 = plt.subplots()
set_axes(ax2,plt,r'Temperature ($K$)',r'$C_{l} (nm)$') 

for pdi_val in pdi_arr:
    if pdi_val == 'expts':
        pdileg = 'PDI: Experimental Distribution'
    else:
        pdileg = 'PDI: ' + str(pdi_val)
    fname = '/Rgscaling_'+str(pdi_val)+'.dat'
    if not os.path.exists(anaout_dir + '/' + fname):
        print('ERR: '+fname+' does not exist in ' + anaout_dir)
        continue
    
    df=pd.read_table(anaout_dir + '/' + fname)
    print('Plotting', pdi_val)
    ax2.scatter(x=df['T'],y=df['nu'],marker=mrk_arr[pdiindx],\
                color=clr_arr[pdiindx])
    ax3.scatter(x=df['T'],y=df['b'],marker=mrk_arr[pdiindx],\
                color=clr_arr[pdiindx])
    
fig2.savefig(figout_dir + '/'+'nu_allpdi.png',dpi=fig2.dpi)
fig2.savefig(figout_dir + '/'+'nu_allpdi.eps',format='eps')
plt.close(fig2)

fig3.savefig(figout_dir + '/'+'perlen_allpdi.png',dpi=fig3.dpi)
fig3.savefig(figout_dir + '/'+'perlen_allpdi.eps',format='eps')
plt.close(fig3)
#----------------------End of Rg2 Analysis------------------------------
