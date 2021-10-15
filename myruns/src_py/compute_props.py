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
import os
import sys
import re
import shutil
import glob
import math
import subprocess
#------------------------------------------------------------------

# Compute Tg
def tg_compute(df,ax):

    tvals = np.array(df['Temp']); den_vals = np.array(df['SV_NPT'])
    lden  = len(tvals); 
    if lden < 7:
        print("ERR:Requires at least 7 temperature points to fit")
        return 0

    clf = LinearRegression(fit_intercept = True) 
    toterr = 99999999999 

    for cnt in range(2,lden-1):

        clf.fit(tvals[0:cnt].reshape(-1,1),den_vals[0:cnt])
        predlo = clf.predict(tvals[0:cnt].reshape(-1,1))
        errlo  = metrics.mean_squared_error(den_vals[0:cnt],predlo)

        clf.fit(tvals[cnt:lden].reshape(-1,1),den_vals[cnt:lden])
        predhi = clf.predict(tvals[cnt:lden].reshape(-1,1))
        errhi  = metrics.mean_squared_error(den_vals[cnt:lden],predhi)
        if errlo + errhi < toterr:
            refcnt = cnt; toterr = errlo + errhi
            refpredlo = predlo; refpredhi = predhi
    ax.plot(tvals[0:refcnt],refpredlo,'--',color='red',linewidth=2)
    ax.plot(tvals[refcnt:lden],refpredhi,'--',color='red',linewidth=2)
    return tvals[refcnt]
#------------------------------------------------------------------ 

# if __name__
if __name__ == '__main__':
    main()
#------------------------------------------------------------------







