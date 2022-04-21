# Analyze and plot Hydrogen Bonds
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
temp_dt   = 50  # Temperature dt
pdi_arr   = [1.0,1.8,3.0,'expts']
mark_arr  = ['o','d','s']
nchains   = 20
anadir_head = 'all_hbs' # change this for different analysis
start_frac  = 0.3 # starting point for averaging
#------------------------------------------------------------------

# Global arrays
denarr = np.arange(temp_min,temp_max,5*temp_dt) #plotting time data
#------------------------------------------------------------------

# Generate output directories
anaout_dir = anaout_dir + '/hb_results' # result_outputs
if not os.path.isdir(anaout_dir):
    os.mkdir(anaout_dir)
#------------------------------------------------------------------

#----------------Begin Analyzing HB data -----------------------------
print("Analyzing HB data ...")
for pdi_val in pdi_arr:

    # Set/Check directories
    simout_dir = all_dir + '/results_pdi_' + str(pdi_val)
    if not os.path.isdir(simout_dir):
        print("ERR: " +  simout_dir + " does not exist")       
        continue
    
    # Output for each case
    fall_out = open(anaout_dir +'/allhbintradata_'+str(pdi_val)+'.dat','w')
    fall_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %('Temp','hb_R1','hb_R2',\
                                                'hb_R3','hb_R4','TotCase'))

    fall2_out = open(anaout_dir +'/allhbinterdata_'+str(pdi_val)+'.dat','w')
    fall2_out.write('%s\t%s\t%s\t%s\t%s\t%s\n' %('Temp','hb_R1','hb_R2',\
                                                 'hb_R3','hb_R4','TotCase'))


    # Averaged outputs
    fc_hb = open(anaout_dir +'/HBdata_'+str(pdi_val)+'.dat','w')
    fc_hb.write('%s\t%s\t%s\t%s\t%s\n' \
                %('Temp','Intra_HB','Inter_HB','Tot_HB','Ncases'))

    pdiflag = -1                # to account for empty directories
    
    # Temperature loops and averaging
    for tval in range(temp_min,temp_max,temp_dt): # loop in temp
        temp_leg  = str(tval)
        pdiflag = 1
        intrach_hb = 0; interch_hb = 0; ncases_pertemp = 0; tot_hb = 0
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
            intra_list_of_files = glob.glob(wdir + '/hb_nptmain_intra_*.xvg')
            if intra_list_of_files == []:
                fall_out.write('%s\t' %('N/A')); fall2_out.write('%s\t' %('N/A'))
                print("Intra HB files do not exist for ", tval)
                continue

            inter_list_of_files = glob.glob(wdir + '/hb_nptmain_inter_*.xvg')
            if inter_list_of_files == []:
                fall2_out.write('%s\t' %('N/A'))
                print("Inter HB files do not exist for ", tval)
                continue

            if len(intra_list_of_files) != nchains or \
               len(inter_list_of_files) != nchains:
                fall_out.write('%s\t' %('DiffNchains'))
                fall2_out.write('%s\t' %('DiffNchains'))
                print('ERR: Mismatch in number of analysis file and input',\
                      nchains, len(inter_list_of_files))
                continue

            fcase_hb = open(anaout_dir + '/hb_casenum_' + str(run_arr[casenum])+ \
                      '_T_' + str(tval) + '.dat','w')
            fcase_hb.write('%s\t%s\t%s\t%s\t%s\n' \
                           %('ChainID','Nmons','IntraHB','InterHB','TotHB'))

            case_intra = 0; case_inter = 0
            if not os.path.exists(wdir + '/chainlist.dat'):
                mon_arr = np.full((nchains),20)

                #fall_out.write('%s\t' %('Nochainlist'))
                #fall2_out.write('%s\t' %('Nochainlist'))
                #fall3_out.write('%s\t' %('Nochainlist'))
                print('ERR: chainlist.dat not found')
                #continue
            else:
                mon_arr = ret_mons(wdir + '/chainlist.dat')

            for fyle in intra_list_of_files: # intra chain loop
                chid = ret_chid(fyle)

                fintra = wdir + '/hb_nptmain_intra_'+str(chid)+ '.xvg'
                finter = wdir + '/hb_nptmain_inter_'+str(chid)+ '.xvg'

                if not os.path.exists(fintra) or not \
                   os.path.exists(finter):
                    print("FATAL ERR: Intra/Inter HB files of same chain not found")
                   
                
                # Open and parse intra HB file
                with open(fintra) as fin:
                    lines = (line.lstrip() for line in fin \
                             if not line.lstrip().startswith('#') and \
                             not line.lstrip().startswith('@'))
                    data  = np.loadtxt(lines) #end with open(fintra)
                    
                l1 = int(start_frac*len(data[:1]))
                l2 = len(data[:,1])
                intra = np.average((data[l1:l2,1]))
                case_intra += intra

                # Open and parse inter file
                with open(finter) as fin:
                    lines = (line.lstrip() for line in fin \
                             if not line.lstrip().startswith('#') and \
                             not line.lstrip().startswith('@'))
                    data  = np.loadtxt(lines) #end with open(finter)
            
                l1 = int(start_frac*len(data[:1]))
                l2 = len(data[:,1])
                inter = np.average((data[l1:l2,1]))
                case_inter += inter

                fcase_hb.write('%g\t%g\t%g\t%g\t%g\n'\
                               %(chid,int(mon_arr[chid-1]),intra,inter,\
                                 intra+inter))
                 
            # end for fyle in inter/intra chain loops

            fcase_hb.close() # close writing for each case 
            intrach_hb += case_intra/nchains
            interch_hb += case_inter/nchains
            tot_hb     += case_intra/nchains + case_inter/nchains
            ncases_pertemp += 1 

            # Write average intra/inter hbonds
            fall_out.write('%g\t' %(case_intra/nchains))
            fall2_out.write('%g\t' %(case_inter/nchains)) 
            # end for casenum in len(range(run_arr))
        
        # Do NOT continue if zero cases are found
        if ncases_pertemp == 0:
            fall_out.write('%g\n' %(ncases_pertemp))
            fall2_out.write('%g\t' %(ncases_pertemp))
            continue

        # Write averages over different cases for each temp
        intrach_hb /= ncases_pertemp
        interch_hb /= ncases_pertemp
        tot_hb     /= ncases_pertemp
        fc_hb.write('%g\t%g\t%g\t%g\t%d\n' %(tval,intrach_hb,\
                                             interch_hb,tot_hb,ncases_pertemp))
        fall_out.write('%g\n' %(ncases_pertemp))
        fall2_out.write('%g\n' %(ncases_pertemp))

    # Close temp files
    fc_hb.close();fall_out.close(); fall2_out.close()

    # Do NOT continue if no PDI-temps are found
    if pdiflag == -1: 
        continue

    # Plot HBintra/inter-temperature data for each PDI
    df=pd.read_table(anaout_dir +'/HBdata_'+str(pdi_val)+'.dat')
    plot_allhb(df,figout_dir,pdi_val) #end of PDI loop
    
#------- Plot HB-Temp data for all PDI values together--------------------
print("Plotting all HB data as a function of PDI..")
fig2, ax2 = plt.subplots()
set_axes(ax2,plt,r'Temperature ($K$)',r'#HB$_{\rm{Intra}}$')

fig3, ax3 = plt.subplots()
set_axes(ax3,plt,r'Temperature ($K$)',r'#HB$_{\rm{Inter}}$')

fig4, ax4 = plt.subplots()
set_axes(ax4,plt,r'Temperature ($K$)',r'#HB$_{\rm{Total}}$')

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
    ax2.scatter(x=df['Temp'],y=df['Intra_HB'],label=pdileg)
    ax3.scatter(x=df['Temp'],y=df['Inter_HB'],label=pdileg)
    ax4.scatter(x=df['Temp'],y=df['Tot_HB'],label=pdileg)

fig2.savefig(figout_dir + '/'+'hbintra_allpdi.png')
fig2.savefig(figout_dir + '/'+'hbintra_allpdi.eps',format='eps')
plt.legend(loc=0)
plt.close(fig2)


fig3.savefig(figout_dir + '/'+'hbinter_allpdi.png')
fig3.savefig(figout_dir + '/'+'hbinter_allpdi.eps',format='eps')
ax3.legend(loc=0)
plt.close(fig3)

fig4.savefig(figout_dir + '/'+'hbtot_allpdi.png')
fig4.savefig(figout_dir + '/'+'hbtot_allpdi.eps',format='eps')
ax4.legend(loc=0)
plt.close(fig4)
#----------------------End of Hbond Analysis------------------------------
