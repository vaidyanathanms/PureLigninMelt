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

#----------------Begin Analyzing HB data -----------------------------
print("Analyzing MSD data ...")
for pdi_val in pdi_arr:

    # Set/Check directories
    simout_dir = all_dir + '/results_pdi_' + str(pdi_val)
    if not os.path.isdir(simout_dir):
        print("ERR: " +  simout_dir + " does not exist")       
        continue
    
    pdiflag = -1                # to account for empty directories
    
    # Temperature loops and averaging
    for tval in range(temp_min,temp_max,temp_dt): # loop in temp
        temp_leg  = str(tval)
        pdiflag = 1;  msdtime = 0; ncases_pertemp = 0

        # Output for each case - different from HBonds since this is
        # at a specific temperatture
        fall_out = open(anaout_dir +'/allmsddata_T_'+str(tval)+\
                        str(pdi_val)+'.dat','w')
        fall_out.write('%s\t%g\t%s\n' %('MSD @t = ', tref, ' ps'))
        fall_out.write('%s\t%s\t%s\t%s\n' %('CaseNum','ChID','NMons','MSD'))

        for casenum in range(len(run_arr)): # loop in runarr
            wdir = simout_dir + '/run_' + str(run_arr[casenum]) +\
                '/T_' + str(tval) + '/' + anadir_head
            if not os.path.isdir(wdir):
                fall_out.write('%s\t' %('N/A'))
                print("ERR: " + wdir + " does not exist")
                continue
            
            print("Analyzing: ", pdi_val,tval,run_arr[casenum])
            
            # Check file(s)
            msd_list_of_files = glob.glob(wdir + '/msd_nptmain_*.xvg')
            if msd_list_of_files == []:
                fall_out.write('%s\t' %('N/A'))
                print("MSD files do not exist for ", tval)
                continue

            if len(msd_list_of_files) != nchains:
                fall_out.write('%s\t' %('DiffNchains'))
                print('ERR: Mismatch in number of analysis file and input',\
                      nchains, len(msd_list_of_files))
                continue

            fcase_msd = open(anaout_dir + '/msd_pdi_' + str(pdi_val) + '_casenum_'\
                            + str(run_arr[casenum])+ '_T_' + str(tval)
                            + '.dat','w')
            fcase_msd.write('%s\t%g\t%s\n' %('MSD @t = ', tref, ' ps'))
            fcase_msd.write('%s\t%s\t%s\n' %('ChainID','Nmons','MSD'))

            case_msd = 0
            if not os.path.exists(wdir + '/chainlist.dat'):
                fall_out.write('%s\t' %('Nochainlist'))
                print('ERR: chainlist.dat not found')
                continue
            else:
                mon_arr = ret_mons(wdir + '/chainlist.dat')

            for fyle in msd_list_of_files: # msd chain loop
                chid = ret_chid(fyle) 
                fmsd = wdir + '/msd_nptmain_'+str(chid)+ '.xvg' #inp file

                if not os.path.exists(fmsd):
                    print("FATAL ERR: MSD files of chain " + str(chid)\
                          + " not found")
                    fmsd.close()
                    continue
            
                # Open and parse msd file
                with open(fmsd) as fin:
                    lines = (line.lstrip() for line in fin \
                             if not line.lstrip().startswith('#') and \
                             not line.lstrip().startswith('@'))
                    data  = np.loadtxt(lines) #end with open(fmsd)
                    
                msd = cprop.msd_tref(data,tref)                   
                case_msd += msd
                fcase_msd.write('%g\t%g\t%g\n' %(chid,mon_arr[chid],case_msd))
                fall_out.write('%g\t%g\t%g\t%g\n' %(casenum,chid,\
                                                        mon_arr[chid],case_msd))
                 
            # end for fyle in inter/msd chain loops
            fcase_msd.close() # close writing for each case 

        
        # Do NOT continue if zero cases are found
        if ncases_pertemp == 0:
            fall_out.close()
            continue

    # Close temp files
    fall_out.close()

    # Do NOT continue if no PDI-temps are found
    if pdiflag == -1: 
        continue
#------- Plot HB-Temp data for all PDI values together--------------------
print("Plotting all HB data as a function of PDI..")
fig2, ax2 = plt.subplots()
set_axes(ax2,plt,r'Temperature ($K$)',r'#HB$_{\rm{Msd}}$')

fig3, ax3 = plt.subplots()
set_axes(ax3,plt,r'Temperature ($K$)',r'#HB$_{\rm{Inter}}$')

fig4, ax4 = plt.subplots()
set_axes(ax4,plt,r'Temperature ($K$)',r'#HB$_{\rm{Total}}$')

ymaxmsd = 0; ymaxinter=0; ymaxtotal=0
yminmsd = 1000; ymininter=1000; ymintotal=1000
for pdi_val in pdi_arr:
    if pdi_val == 'expts':
        pdileg = 'PDI: Experimental Distribution'
    else:
        pdileg = 'PDI: ' + str(pdi_val)
    fname = '/HBdata_'+str(pdi_val)+'.dat'
    if not os.path.exists(anaout_dir + '/' + fname):
        print('ERR: '+fname+' does not exist in ' + anaout_dir)
        continue
    
    df=pd.read_table(anaout_dir + '/' + fname)
    print('Plotting', pdi_val)
    ax2.scatter(x=df['Temp'],y=df['Msd_HB'],label=pdileg)
    ax3.scatter(x=df['Temp'],y=df['Inter_HB'],label=pdileg)
    ax4.scatter(x=df['Temp'],y=df['Tot_HB'],label=pdileg)
    ymaxmsd, yminmsd = axlims(ymaxmsd,df['Msd_HB'].max(),\
                                  yminmsd,df['Msd_HB'].min()) 
    ymaxinter, ymininter = axlims(ymaxinter,df['Inter_HB'].max(),\
                                  ymininter,df['Inter_HB'].min()) 
    ymaxtotal, ymintotal = axlims(ymaxtotal,df['Tot_HB'].max(),\
                                  ymintotal,df['Tot_HB'].min()) 


fig2.savefig(figout_dir + '/'+'hbmsd_allpdi.png')
fig2.savefig(figout_dir + '/'+'hbmsd_allpdi.eps',format='eps')
plt.legend(loc='upper right')
ax2.set_ylim([0.9*yminmsd, 1.2*ymaxmsd])
plt.close(fig2)

fig3.savefig(figout_dir + '/'+'hbinter_allpdi.png')
fig3.savefig(figout_dir + '/'+'hbinter_allpdi.eps',format='eps')
ax2.set_ylim([0.9*ymininter, 1.2*ymaxinter])
plt.legend(loc='upper right')
plt.close(fig3)

fig4.savefig(figout_dir + '/'+'hbtot_allpdi.png')
fig4.savefig(figout_dir + '/'+'hbtot_allpdi.eps',format='eps')
ax4.set_ylim([0.9*ymintotal, 1.2*ymaxtotal])
plt.legend(loc='upper right')
plt.close(fig4)
#----------------------End of Hbond Analysis------------------------------
