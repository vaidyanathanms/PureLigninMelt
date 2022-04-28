# Just plotting all data to avoid redoing the calculations
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
from scipy.stats import sem
#------------------------------------------------------------------

#Input data
pdi_arr   = [1.0,1.8,3.0,3.7,'expts']
temp_arr  = range(300,501,50)

#Input flags
tg_plot    = 0 # Plotting Tg
msd_plot   = 1 # Plotting MSD
rg_plot    = 0 # Plotting Rg
segrg_plot = 0 # Plotting segmental rg
hb_plot    = 0 # Plotting hydrogen bond data

#------- Plot SV-Temp data for all PDI values together ------------
while tg_plot: # Stupid Python won't let me break if loops easily
    
    print("Plotting Tg data")
    anaout_dir = anaout_dir + '/dens_results' # result_outputs
    if not os.path.isdir(anaout_dir):
        print("ERROR: ", + anaout_dir + " not found")
        break

    fig, ax = plt.subplots()
    set_axes(ax,plt,r'Temperature ($K$)',\
             r'Specific Volume ($cm^{3}/g$)')
    ymaxref = 0; yminref = 1000; indx=0
    
    for pdi_val in pdi_arr:
        if pdi_val == 'expts':
            pdileg = 'PDI: Experimental Distribution'
        else:
            pdileg = 'PDI: ' + str(pdi_val)
        fplot = '/tgdata_'+str(pdi_val)+'.dat'
        if not os.path.exists(anaout_dir + '/' + fplot):
            print('ERR: '+fplot+' does not exist in ' + anaout_dir)
            continue

        df=pd.read_table(anaout_dir + '/' + fplot)
        print('Plotting', pdi_val)
        plt.scatter(x=df['Temp'],y=df['SV_NPT'],marker=mrk_arr[indx]\
                    ,label=pdileg)
        yminref, ymaxref = axlims(yminref,df['SV_NPT'].min(),\
                                  ymaxref,df['SV_NPT'].max()) 
        indx += 1

    plt.legend(loc='upper right')
    ax.set_ylim([0.98*yminref, 1.1*ymaxref])
    fig.savefig(figout_dir + '/'+'svt_allpdi.png',dpi=fig.dpi)
    fig.savefig(figout_dir + '/'+'svt_allpdi.eps',format='eps')
    plt.close(fig)
    tg_plot = 0 
#--- End plotting Tg data -----------------------------------------

#------- Plot MSD data---------------------------------------------
while msd_plot:

    print("Plotting MSD data")
    anaout_dir = anaout_dir + '/msd_results' # result_outputs
    if not os.path.isdir(anaout_dir):
        print("ERROR: ", + anaout_dir + " not found")
        break

    for tval in temp_arr:

        fig, ax = plt.subplots()
        set_axes(ax,plt,r'$N$',r'MSD ($\AA^{2}$)')
        ymaxref = 0; yminref = 1000; indx=0
        
        for pdi_val in pdi_arr:
            
            if pdi_val == 'expts':
                pdileg = 'PDI: Experimental Distribution'
            else:
                pdileg = 'PDI: ' + str(pdi_val)
            fplot = '/allmsddata_T_'+ str(tval) + '_pdi_' + \
                str(pdi_val)+'.dat'
            
            if not os.path.exists(anaout_dir + '/' + fplot):
                print('ERR: '+fplot+' does not exist in ' + anaout_dir)
                continue

            # Averaging across cases
            print('Averaging', pdi_val, tval)
            df=pd.read_table(anaout_dir + '/' + fplot,skiprows=1)
            mondata = df['NMons']
            msddata = df['MSD']

            favg = open(anaout_dir + '/msdavg_T_'+ str(tval) +\
                        '_pdi_'+str(pdi_val)+'.dat','w')
            favg.write('%s\t%s\t%s\n' %('NMons','AvgMSD','Ncnts'))
            
            av_monarr = []; av_msdarr = []; av_cntarr = []
            for monindx in range(len(mondata)):
                monval = mondata[monindx]; msd = msddata[monindx]
                if math.isnan(monval):
                    continue
                elif monval in av_monarr:
                    av_msdarr[list(av_monarr).index(monval)] += msd
                    av_cntarr[list(av_monarr).index(monval)] += 1
                else:
                    av_monarr = np.append(av_monarr,monval)
                    av_msdarr = np.append(av_msdarr,msd)
                    av_cntarr = np.append(av_cntarr,1)

            for avval in range(len(av_monarr)):
                favg.write('%g\t%g\t%g\n' %(av_monarr[avval],\
                                            av_msdarr[avval],\
                                            av_cntarr[avval]))
            favg.close()

            if len(av_msdarr) != len(av_monarr) or \
               len(av_msdarr) != len(av_cntarr):
                print('FATAL ERR: Sizes not same for pdi/T: ',\
                      pdi_val, tval)
                continue

            #Plotting data
            print('Plotting', pdi_val, tval)
            plt.scatter(av_monarr,100*av_msdarr/av_cntarr,\
                        marker=mrk_arr[indx],label=pdileg) # in Ang
            yminref, ymaxref = axlims(yminref,(av_msdarr/av_cntarr).min(),\
                                      ymaxref,(av_msdarr/av_cntarr).max()) 
            indx += 1
        
#        print(yminref,ymaxref)
        plt.legend(loc='upper right')
        ax.set_ylim([0.9*yminref, 1.2*ymaxref])
        ax.set_yscale('log'); ax.set_xscale('log')
        fig.savefig(figout_dir + '/'+'MSDdist_T_'+str(tval)+ \
                    '.png',dpi=fig.dpi)
        fig.savefig(figout_dir + '/'+'MSDdist_T_'+str(tval)+ \
                    '.eps',format='eps')
        plt.close(fig)
    msd_plot = 0
#----------------------End MSD plots--------------------------------------

#------- Plot HB-Temp data for all PDI values together--------------------
while hb_plot:

    print("Plotting all HB data")
    anaout_dir = anaout_dir + '/msd_results' # result_outputs
    if not os.path.isdir(anaout_dir):
        print("ERROR: ", + anaout_dir + " not found")
        break
    
    fig2, ax2 = plt.subplots()
    set_axes(ax2,plt,r'Temperature ($K$)',r'#HB$_{\rm{Intra}}$')

    fig3, ax3 = plt.subplots()
    set_axes(ax3,plt,r'Temperature ($K$)',r'#HB$_{\rm{Inter}}$')

    fig4, ax4 = plt.subplots()
    set_axes(ax4,plt,r'Temperature ($K$)',r'#HB$_{\rm{Total}}$')

    ymaxintra = 0; ymaxinter=0; ymaxtotal=0
    yminintra = 1000; ymininter=1000; ymintotal=1000
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
        ymaxintra, yminintra = axlims(ymaxintra,df['Intra_HB'].max(),\
                                      yminintra,df['Intra_HB'].min()) 
        ymaxinter, ymininter = axlims(ymaxinter,df['Inter_HB'].max(),\
                                      ymininter,df['Inter_HB'].min()) 
        ymaxtotal, ymintotal = axlims(ymaxtotal,df['Tot_HB'].max(),\
                                      ymintotal,df['Tot_HB'].min()) 


    plt.legend(loc='upper right')
    ax2.set_ylim([0.9*yminintra, 1.2*ymaxintra])
    fig2.savefig(figout_dir + '/'+'hbintra_allpdi.png')
    fig2.savefig(figout_dir + '/'+'hbintra_allpdi.eps',format='eps')
    plt.close(fig2)

    ax3.set_ylim([0.9*ymininter, 1.2*ymaxinter])
    plt.legend(loc='upper right')
    fig3.savefig(figout_dir + '/'+'hbinter_allpdi.png')
    fig3.savefig(figout_dir + '/'+'hbinter_allpdi.eps',format='eps')
    plt.close(fig3)

    ax4.set_ylim([0.9*ymintotal, 1.2*ymaxtotal])
    plt.legend(loc='upper right')
    fig4.savefig(figout_dir + '/'+'hbtot_allpdi.png')
    fig4.savefig(figout_dir + '/'+'hbtot_allpdi.eps',format='eps')
    plt.close(fig4)
    hb_plot = 0
#----------------------End Hbond plots-------------------------------------
