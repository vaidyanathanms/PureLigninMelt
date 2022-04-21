# All path and color inputs for analyzing/plotting data

#Import Modules
import os
import matplotlib.pyplot as plt

# Directory paths
main_dir = os.getcwd() # current dir
all_dir  = '../../../sim_results' # analysis dir
if not os.path.isdir(all_dir):
    print("FATAL ERROR: ", all_dir, " not found")
    exit("Check scratch directory path")
anaout_dir = '../../../analyzed_results'
if not os.path.isdir(anaout_dir):
    os.mkdir(anaout_dir)
figout_dir = '../../../figures'
if not os.path.isdir(figout_dir):
    os.mkdir(figout_dir)  
#------------------------------------------------------------------

# Color/line data; figure defaults
orange = '#FFA500'; dark_g = '#006400'; brown = '#8B4513'
clr_arr = [dark_g,orange,brown,'b','k','m']
mrk_arr = ['o','d','s']
lne_arr = ['-','--']
plt.rcParams['figure.dpi'] = 1200
plt.rc('legend',fontsize=16) # fontsize of the legends
plt.rcParams.update({'font.size': 16}) # other fonts
plt.rcParams.update({'figure.autolayout': True})
#plt.rcParams['font.family'] = 'Helvetica'
#plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
#plt.rcParams.update({'font.family': 'Times New Roman'})
#------------------------------------------------------------------

# Input path data
inp_type  = 'melts' # melts, solvents, cosolvents
disperse  = 'poly' # mono/poly; only for melts
biom_arr  = ['WT'] # biomass type arr
otyp_arr  = ['None']  # solvent arr for solvents/cosolvents
otyp_leg  = ['None']  # legend for solvents/cosolvents
solv_type = 'None'
#------------------------------------------------------------------
