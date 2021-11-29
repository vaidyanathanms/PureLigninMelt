# To compute Tg

# Import modules
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import os
import sys
import re
import shutil
import glob
import math
import subprocess
from plt_aux import *
from compute_props import *
#--------------------------------------------------------------------------------

fig_dir = 'path_to_output_directory'
tg_fyl = 'input_file.dat'
pdi_val = 1.0 #use this as a dummy
df=pd.read_table(tg_fyl) # read the file
figa, axa = plot_tg(df,fig_dir,pdi_val) #figure settings
tgval = compute_tg(df,axa) #compute tgval and plot
