# Plot all data
# Version: Apr-28-2022
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
temp_arr  = range(300,501,10) # for temperature specific plots

#Input flags
tg_plot    = 1 # Plotting Tg
msd_plot   = 0 # Plotting MSD
rg_plot    = 0 # Plotting average Rg
bnu_plot   = 0 # Plotting segmental Rg
rgNM_plot  = 0 # Plotting segmental Rg as a function of deg of poly
hb_plot    = 0 # Plotting hydrogen bond data
shape_plot = 0 # Plotting shape factor

#--------Plot SV-Temp data for all PDI values together ------------
while tg_plot: # Stupid Python won't let me break if loops easily
    
    print("Plotting Tg data")
    anaout1_dir = anaout_dir + '/dens_results' # result_outputs
    if not os.path.isdir(anaout1_dir):
        print("ERROR: ", + anaout1_dir + " not found")
        break

    # Plot all dataa
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
        if not os.path.exists(anaout1_dir + '/' + fplot):
            print('ERR: '+fplot+' does not exist in ' + anaout1_dir)
            continue
        
        print('Plotting', pdi_val)
        df=pd.read_table(anaout1_dir + '/' + fplot)
        plt.errorbar(x=df['Temp'],y=df['SV_NPT'],yerr=df['SV_err'],\
                     marker=mrk_arr[indx],label=pdileg,\
                     capsize=5,linestyle='None')
        yminref, ymaxref = axlims(yminref,df['SV_NPT'].min(),\
                                  ymaxref,df['SV_NPT'].max()) 
        indx += 1

    plt.legend(loc='upper left')
#    print(ymaxref)
    ax.set_ylim([0.98*yminref, 1.05*ymaxref])
    fig.savefig(figout_dir + '/'+'svt_allpdi.png',dpi=fig.dpi)
    fig.savefig(figout_dir + '/'+'svt_allpdi.eps',format='eps')
    plt.close(fig)
    tg_plot = 0     

    # # Plot average Tg
    # fig2, ax2 = plt.subplots()
    # set_axes(ax2,plt,r'Temperature ($K$)',\
    #          r'Specific Volume ($cm^{3}/g$)')
    # ymaxref = 0; yminref = 1000; indx=0
    
    # for pdi_val in pdi_arr:
    #     if pdi_val == 'expts':
    #         pdileg = 'PDI: Experimental Distribution'
    #     else:
    #         pdileg = 'PDI: ' + str(pdi_val)
    #     fplot = '/tgdata_'+str(pdi_val)+'.dat'
    #     if not os.path.exists(anaout1_dir + '/' + fplot):
    #         print('ERR: '+fplot+' does not exist in ' + anaout1_dir)
    #         continue
        
    #     print('Plotting', pdi_val)
    #     df=pd.read_table(anaout1_dir + '/' + fplot)
    #     plt.scatter(x=df['Temp'],y=df['SV_NPT'],marker=mrk_arr[indx]\
    #                 ,label=pdileg)
    #     yminref, ymaxref = axlims(yminref,df['SV_NPT'].min(),\
    #                               ymaxref,df['SV_NPT'].max()) 
    #     indx += 1

    # plt.legend(loc='upper right')
    # ax.set_ylim([0.95*yminref, 1.2*ymaxref])
    # fig.savefig(figout_dir + '/'+'svt_allpdi.png',dpi=fig.dpi)
    # fig.savefig(figout_dir + '/'+'svt_allpdi.eps',format='eps')
    # plt.close(fig)
    


#--- End plotting Tg data -----------------------------------------

#--------Plot MSD data for all PDI together------------------------
while msd_plot:

    print("Plotting MSD data")
    anaout1_dir = anaout_dir + '/msd_results' # result_outputs
    if not os.path.isdir(anaout1_dir):
        print("ERROR: ", + anaout1_dir + " not found")
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
            
            if not os.path.exists(anaout1_dir + '/' + fplot):
                print('ERR: '+fplot+' does not exist in ' + anaout1_dir)
                continue

            # Averaging across cases
            print('Averaging', pdi_val, tval)
            df=pd.read_table(anaout1_dir + '/' + fplot,skiprows=1)
            mondata = df['NMons']
            msddata = df['MSD']
            favg = open(anaout1_dir + '/msdavg_T_'+ str(tval) +\
                        '_pdi_'+str(pdi_val)+'.dat','w')
            favg.write('%s\t%s\t%s\t%s\n' %('NMons','AvgMSD',\
                                            'ErrMSD','Ncnts'))
            
            print('Plotting', pdi_val, tval)
            monplot,msdplot,errplot,cntall = comp_bin_ave_sem(mondata,msddata)
            plt.errorbar(x=monplot,y=100*msdplot,yerr=100*errplot,\
                         marker=mrk_arr[indx],label=pdileg,\
                         capsize=5,linestyle='None') #100 for nm to A
            yminref, ymaxref = axlims(yminref,(msdplot-errplot).min(),\
                                      ymaxref,(msdplot+errplot).max()) 
            indx += 1

        plt.legend(fontsize=10,loc='upper right')
#        plt.legend(bbox_to_anchor=(0.45,0.65),fontsize=10)
        ax.set_ylim([95*yminref, 120*ymaxref])
        ax.set_yscale('log'); ax.set_xscale('log')
        fig.savefig(figout_dir + '/'+'MSDdist_T_'+str(tval)+ \
                    '.png',dpi=fig.dpi)
        fig.savefig(figout_dir + '/'+'MSDdist_T_'+str(tval)+ \
                    '.eps',format='eps')
        plt.close(fig)
            
    msd_plot = 0
#----------------------End MSD plots--------------------------------------

#--------Plot Rg2-Temp data for all PDI values together-------------------
while rg_plot:
    print("Plotting Rg data")
    anaout1_dir = anaout_dir + '/rg_results' # result_outputs
    if not os.path.isdir(anaout1_dir):
        print("ERROR: ", + anaout1_dir + " not found")
        break

    fig, ax = plt.subplots()
    set_axes(ax,plt,r'Temperature ($K$)',\
             r'$\langle R_{g}^{2} \rangle$ ($nm^{2}$)')
    ymaxref = 0; yminref = 1000; indx=0
        
    for pdi_val in pdi_arr:
        if pdi_val == 'expts':
            pdileg = 'PDI: Experimental Distribution'
        else:
            pdileg = 'PDI: ' + str(pdi_val)
        fplot = '/Rgdata_'+str(pdi_val)+'.dat'
        if not os.path.exists(anaout1_dir + '/' + fplot):
            print('ERR: '+fplot+' does not exist in ' + anaout1_dir)
            continue

        print('Plotting', pdi_val)
        df=pd.read_table(anaout1_dir + '/' + fplot)
        plt.scatter(x=df['Temp'],y=df['<Rg^2>'],marker=mrk_arr[indx],\
                    label=pdileg); indx+=1
        yminref, ymaxref = axlims(yminref,df['<Rg^2>'].min(),\
                                  ymaxref,df['<Rg^2>'].max())
    
    plt.legend(loc='upper right')
    ax.set_ylim([0.95*yminref, 1.2*ymaxref])
    fig.savefig(figout_dir + '/'+'rg2_allpdi.png',dpi=fig.dpi)
    fig.savefig(figout_dir + '/'+'rg2_allpdi.eps',format='eps')
    plt.close(fig)
    rg_plot = 0
#----------------------End Rg2 plots--------------------------------------

#--------Plot b/nu-Temp data for all PDI values together------------------
while bnu_plot:

    print("Plotting Rg Scaling Coefficients")
    anaout1_dir = anaout_dir + '/segrg_results' # result_outputs
    if not os.path.isdir(anaout1_dir):
        print("ERROR: ", + anaout_dir + " not found")
        break

    fig, ax = plt.subplots()
    set_axes(ax,plt,r'Temperature ($K$)',r'$\nu$') 
    ymaxnuref = 0; yminnuref = 1000
    
    fig2, ax2 = plt.subplots()
    set_axes(ax2,plt,r'Temperature ($K$)',r'$C_{l}$ ($\AA$)') 
    ymaxbref = 0; yminbref = 1000; indx = 0
   
    for pdi_val in pdi_arr:
        if pdi_val == 'expts':
            pdileg = 'PDI: Experimental Distribution'
        else:
            pdileg = 'PDI: ' + str(pdi_val)
        fplot = '/Rgscaling_'+str(pdi_val)+'.dat'
        if not os.path.exists(anaout1_dir + '/' + fplot):
            print('ERR: '+fplot+' does not exist in ' + anaout1_dir)
            continue

        print('Plotting', pdi_val)
        df=pd.read_table(anaout1_dir + '/' + fplot)
        ax.scatter(x=df['Temp'],y=df['nu'],marker=mrk_arr[indx]\
                   ,label = pdileg); indx += 1
        yminnuref, ymaxnuref = axlims(yminnuref,df['nu'].min(),\
                                  ymaxnuref,df['nu'].max())
        ax2.scatter(x=df['Temp'],y=df['b'],label = pdileg)
        yminbref, ymaxbref = axlims(yminbref,df['b'].min(),\
                                  ymaxbref,df['b'].max())

        
    ax.legend(loc='upper right')
    ax.set_ylim([0.95*yminnuref, 1.2*ymaxnuref])
    fig.savefig(figout_dir + '/'+'nu_allpdi.png',dpi=fig.dpi)
    fig.savefig(figout_dir + '/'+'nu_allpdi.eps',format='eps')
    plt.close(fig)

    ax2.legend(loc='upper right')
    ax2.set_ylim([0.95*yminbref, 1.2*ymaxbref])
    fig2.savefig(figout_dir + '/'+'perlen_allpdi.png',dpi=fig2.dpi)
    fig2.savefig(figout_dir + '/'+'perlen_allpdi.eps',format='eps')
    plt.close(fig2)
    bnu_plot = 0
#----------------------End scaling plots----------------------------------

#--------Plot Rg as a function of degree of polymerization ---------------
while rgNM_plot:

    print("Plotting Rg as a function of degree of polymerization")
    anaout1_dir = anaout_dir + '/segrg_results' # result_outputs
    if not os.path.isdir(anaout1_dir):
        print("ERROR: ", + anaout_dir + " not found")
        break
    Nmax = 25 #maximum N at which the plot is cut if lxy>length of polymer
            
    for tempval in temp_arr:

        fig, ax = plt.subplots()
        set_axes(ax,plt,r'$N$',r'$R_g$ ($\AA$)') 
        ymaxref = 0; yminref = 100000; indx =0
    
        for pdi_val in pdi_arr:
            if pdi_val == 'expts':
                pdileg = 'PDI: Experimental Distribution'
            else:
                pdileg = 'PDI: ' + str(pdi_val)
            fplot = '/allRgdata_T_'+str(tempval)+'_PDI_'+\
                str(pdi_val)+'.dat'
            if not os.path.exists(anaout1_dir + '/' + fplot):
                print('ERR: '+fplot+' does not exist in ' + anaout1_dir)
                continue

            print('Plotting T/PDI: ', tempval, pdi_val)
            df=pd.read_table(anaout1_dir + '/' + fplot)
            xo = df['N']; yo = df['RgMean']; yoerr = df['RgSEM']
            lxy = min(len(xo),Nmax)
            plt.errorbar(x=xo[:lxy],y=yo[:lxy],yerr=yoerr[:lxy],\
                         marker=mrk_arr[indx],label=pdileg,\
                         capsize=5,linestyle='None')
            yminref, ymaxref = axlims(yminref,yo[:lxy].min(),\
                                      ymaxref,yo[:lxy].max())
            indx += 1
        ax.set_yscale('log'); ax.set_xscale('log')
        ax.legend(loc='upper left')
        ax.set_ylim([0.95*yminref, 1.2*ymaxref])
        fig.savefig(figout_dir + '/'+'RgN_T_' + str(tempval)\
                    +'.png',dpi=fig.dpi)
        fig.savefig(figout_dir + '/'+'RgN_T_' + str(tempval)\
                    +'.eps',format='eps')
        plt.close(fig)
    rgNM_plot = 0
#----------------------End scaling plots----------------------------------

#-------Plot Kappa-Temp data for all PDI values together------------------
while shape_plot:

    print("Plotting shape data")
    anaout1_dir = anaout_dir + '/shape_results' # result_outputs
    if not os.path.isdir(anaout1_dir):
        print("ERROR: ", + anaout1_dir + " not found")
        break

    for tval in temp_arr:
        
        fig, ax = plt.subplots()
        set_axes(ax,plt,r'$N$',r'$\kappa$')
        ymaxref = 0; yminref = 1000; indx=0
    
        for pdi_val in pdi_arr:
            if pdi_val == 'expts':
                pdileg = 'PDI: Experimental Distribution'
            else:
                pdileg = 'PDI: ' + str(pdi_val)

            fplot = '/shape_T_'+ str(tval) + '_pdi_' + \
                str(pdi_val)+'.dat'
            if not os.path.exists(anaout1_dir + '/' + fplot):
                print('ERR: '+fplot+' does not exist in ' + anaout1_dir)
                continue

            # Averaging across cases
            print('Averaging', pdi_val, tval)
            df=pd.read_table(anaout1_dir + '/' + fplot)
            mondata = df['NMons']
            kapdata = df['kappa']

            favg = open(anaout1_dir + '/shapeavg_T_'+ str(tval) +\
                        '_pdi_'+str(pdi_val)+'.dat','w')
            favg.write('%s\t%s\t%s\n' %('NMons','AvgKappa','Ncnts'))
            print('Plotting', pdi_val, tval)
            monplot,shapeplot,errplot,cntall = comp_bin_ave_sem(mondata,kapdata)
            plt.errorbar(x=monplot,y=shapeplot,yerr=errplot,\
                         marker=mrk_arr[indx],label=pdileg,\
                         capsize=5,linestyle='None') #100 for nm to A
            yminref, ymaxref = axlims(yminref,(shapeplot-errplot).min(),\
                                      ymaxref,(shapeplot+errplot).max()) 
            indx += 1

        ax.legend(loc='upper right')
        ax.set_ylim([0.95*yminref, 1.2*ymaxref])
        fig.savefig(figout_dir + '/'+'kapdist_T_' + str(tval) +\
                    '.png',dpi=fig.dpi)
        fig.savefig(figout_dir + '/'+'kapdist_T_' + str(tval) +\
                    '.eps',format='eps')
        plt.close(fig)
    shape_plot = 0
#----------------------End shape plots----------------------------------

#--------Plot HB-Temp data for all PDI values together--------------------
while hb_plot:

    print("Plotting all HB data")
    anaout1_dir = anaout_dir + '/msd_results' # result_outputs
    if not os.path.isdir(anaout1_dir):
        print("ERROR: ", + anaout1_dir + " not found")
        break
    
    fig, ax = plt.subplots()
    set_axes(ax,plt,r'Temperature ($K$)',r'#HB$_{\rm{Intra}}$')

    fig2, ax2 = plt.subplots()
    set_axes(ax2,plt,r'Temperature ($K$)',r'#HB$_{\rm{Inter}}$')

    fig3, ax3 = plt.subplots()
    set_axes(ax3,plt,r'Temperature ($K$)',r'#HB$_{\rm{Total}}$')

    ymaxintra = 0; ymaxinter=0; ymaxtotal=0
    yminintra = 1000; ymininter=1000; ymintotal=1000; indx = 0
    for pdi_val in pdi_arr:
        if pdi_val == 'expts':
            pdileg = 'PDI: Experimental Distribution'
        else:
            pdileg = 'PDI: ' + str(pdi_val)
        fplot = '/HBdata_'+str(pdi_val)+'.dat'
        if not os.path.exists(anaout1_dir + '/' + fplot):
            print('ERR: '+fplot+' does not exist in ' + anaout1_dir)
            continue

        print('Plotting', pdi_val)
        df=pd.read_table(anaout1_dir + '/' + fplot)
        ax.scatter(x=df['Temp'],y=df['Intra_HB'],marker=mrk_arr[indx],\
                   label=pdileg)
        ax2.scatter(x=df['Temp'],y=df['Inter_HB'],marker=mrk_arr[indx],\
                    label=pdileg)
        ax3.scatter(x=df['Temp'],y=df['Tot_HB'],marker=mrk_arr[indx],\
                    label=pdileg); indx += 1
        ymaxintra, yminintra = axlims(ymaxintra,df['Intra_HB'].max(),\
                                      yminintra,df['Intra_HB'].min()) 
        ymaxinter, ymininter = axlims(ymaxinter,df['Inter_HB'].max(),\
                                      ymininter,df['Inter_HB'].min()) 
        ymaxtotal, ymintotal = axlims(ymaxtotal,df['Tot_HB'].max(),\
                                      ymintotal,df['Tot_HB'].min()) 



    ax.set_ylim([0.95*yminintra, 1.2*ymaxintra])
    ax.legend(loc='upper right')
    fig.savefig(figout_dir + '/'+'hbintra_allpdi.png',dpi=fig.dpi)
    fig.savefig(figout_dir + '/'+'hbintra_allpdi.eps',format='eps')
    plt.close(fig)

    ax2.set_ylim([0.95*ymininter, 1.2*ymaxinter])
    ax2.legend(loc='upper right')
    fig2.savefig(figout_dir + '/'+'hbinter_allpdi.png',dpi=fig2.dpi)
    fig2.savefig(figout_dir + '/'+'hbinter_allpdi.eps',format='eps')
    plt.close(fig2)

    ax3.set_ylim([0.95*ymintotal, 1.2*ymaxtotal])
    ax3.legend(loc='upper right')
    fig3.savefig(figout_dir + '/'+'hbtot_allpdi.png',dpi=fig3.dpi)
    fig3.savefig(figout_dir + '/'+'hbtot_allpdi.eps',format='eps')
    plt.close(fig3)
    hb_plot = 0
#----------------------End Hbond plots-------------------------------------
