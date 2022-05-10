# To compute Rg scaling
# To compute Tg from densities
# Use plotdata.py for calling the file
# V1: Oct-12-2021
#-------------------------------------------------------------------

#Import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.linear_model import LinearRegression
from sklearn import metrics
from scipy.optimize import curve_fit
import os
import sys
import re
import shutil
import glob
import math
import subprocess
#------------------------------------------------------------------

# Compute Tg
def compute_tg(df,ax):

    tvals = np.array(df['Temp']); den_vals = np.array(df['SV_NPT'])
    wvals = np.array(df['SV_err']); wvals = 1/wvals

    lden  = len(tvals); 
    if lden < 7:
        print("ERR:Requires at least 7 temperature points to fit")
        return 0

    clf = LinearRegression(fit_intercept = True) 
    toterr = 99999999999 

    for cnt in range(2,lden-1):

        low_fit = clf.fit(tvals[0:cnt].reshape(-1,1),den_vals[0:cnt],\
                          sample_weight = wvals[0:cnt])
        predlo = clf.predict(tvals[0:cnt].reshape(-1,1))
        errlo = np.sum(wvals[0:cnt]*(tvals[0:cnt]-predlo)**2)
#        errlo  = wvals[0:cnt]*metrics.mean_squared_error(den_vals[0:cnt],predlo)
#        print(cnt,errlo)
        
        hi_fit = clf.fit(tvals[cnt:lden].reshape(-1,1),den_vals[cnt:lden],\
                         sample_weight = wvals[cnt:lden])
        predhi = clf.predict(tvals[cnt:lden].reshape(-1,1))
        errhi = np.sum(wvals[cnt:lden]*(tvals[cnt:lden]-predhi)**2)
#        errhi  = wvals[cnt:lden]*metrics.mean_squared_error(den_vals[cnt:lden],predhi)
#        print(cnt,errhi)
        
        if errlo + errhi < toterr:
            refcnt = cnt; toterr = errlo + errhi
            refpredlo = low_fit.predict(tvals[0:cnt+1].reshape(-1,1))
            refpredhi = hi_fit.predict(tvals[cnt-1:lden].reshape(-1,1))
            
    ax.plot(tvals[0:refcnt+1],refpredlo,'--',color='r',\
            linewidth=2,label='Fit_Low')
    ax.plot(tvals[refcnt-1:lden],refpredhi,'--',color='g',\
            linewidth=2,label='Fit_High')
    return tvals[refcnt],low_fit.intercept_,low_fit.coef_,errlo,\
        hi_fit.intercept_,hi_fit.coef_,errhi
#------------------------------------------------------------------ 

# Compute Rgscaling
def compute_rgscaling(df):
    N_all = np.array(df['N']); Rg_all = np.array(df['Rgmean'])
    lenarr = min(math.floor(0.5*len(N_all)),50)
    if lenarr == 0:
        print("ERR: No data found")
        return [0,0], 0
    xarr = N_all[0:lenarr]; yarr = Rg_all[0:lenarr]
    pars, cov = curve_fit(func_powerlaw,xarr,yarr,p0=np.asarray([2, 0.5]))
    stdevs = np.sqrt(np.diag(cov))
    return pars, cov
#------------------------------------------------------------------ 

# Return closest MSD at given tref
def msd_tref(data,reft):
    if data[len(data[:,0])-2,0] < reft:
        print('Simulation time less than reference time')
        print(data[len(data[:,0])-2,0],reft)
        return -1
    else:
        for tind in range(len(data[:,0])):
            if data[tind,0] < reft and data[tind+1,0] >= reft:
                sl = (data[tind+1,1]-data[tind,1])/\
                    (data[tind+1,0]-data[tind,0])
                return data[tind,1] + sl*(reft-data[tind,0])
#------------------------------------------------------------------ 
# Define power law
def func_powerlaw(x,a,b):
    return a*x**b

#------------------------------------------------------------------ 
# if __name__
if __name__ == '__main__':
    main()
#------------------------------------------------------------------







