# All path and color inputs for analyzing/plotting data
# Directory paths
main_dir = os.getcwd() # current dir
scr_dir  = '/lustre/or-scratch/cades-bsd/v0e' # scratch dir
scr_dir  = scr_dir + '/Glassy_lignin'
sh_dir   = '../src_gmx/sh_files'
if not os.path.isdir(scr_dir):
    print("FATAL ERROR: ", scr_dir, " not found")
    exit("Check scratch directory path")
#------------------------------------------------------------------

# Color/line data; figure defaults
orange = '#FFA500'; dark_g = '#006400'; brown = '#8B4513'
clr_arr = [dark_g,orange,brown,'b','k','m']
mrk_arr = ['o','d','s']
lne_arr = ['-','--']
plt.rc('legend',fontsize=16) # fontsize of the legends
plt.rcParams.update({'font.size': 16}) # other fonts
plt.rcParams.update({'figure.autolayout': True})
#plt.rcParams['font.family'] = 'serif'
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
