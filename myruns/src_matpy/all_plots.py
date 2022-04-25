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

#Input flags
tg_plot = 1 # Plotting Tg

#------- Plot SV-Temp data for all PDI values together ------------
while tg_plot:
    
    print("Analyzing Tg data")
    anaout_dir = anaout_dir + '/dens_results' # result_outputs
    if not os.path.isdir(anaout_dir):
        print("ERROR: ", + anaout_dir + " not found")
        break

    fig2, ax2 = plt.subplots()
    set_axes(ax2,plt,r'Temperature ($K$)',r'Specific Volume ($cm^{3}/g$)')
    ymaxref = 0; yminref = 1000
    indx=0
    
    for pdi_val in pdi_arr:
        if pdi_val == 'expts':
            pdileg = 'PDI: Experimental Distribution'
        else:
            pdileg = 'PDI: ' + str(pdi_val)
        fname = '/tgdata_'+str(pdi_val)+'.dat'
        if not os.path.exists(anaout_dir + '/' + fname):
            print('ERR: '+fname+' does not exist in ' + anaout_dir)
            continue

        df=pd.read_table(anaout_dir + '/' + fname)
        print('Plotting', pdi_val)
        plt.scatter(x=df['Temp'],y=df['SV_NPT'],marker=mrk_arr[indx]\
                    ,label=pdileg)
        ymaxref, yminref = axlims(ymaxref,df['SV_NPT'].max(),\
                                  yminref,df['SV_NPT'].min()) 
        indx += 1
        
    plt.legend(loc='upper left')
    ax2.set_ylim([0.9*yminref, 1.2*ymaxref])
    fig2.savefig(figout_dir + '/'+'svt_allpdi.png',dpi=fig2.dpi)
    fig2.savefig(figout_dir + '/'+'svt_allpdi.eps',format='eps')
    plt.close(fig2)
    break # Stupid Python won't let me break if loops easily
#--- End plotting Tg data --------------------------------------------

