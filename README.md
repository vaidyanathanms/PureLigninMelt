# PureLigninMelt
To generate lignin melt systems and compute lignin melt properties.

Author: Vaidyanathan Sethuraman
Email: v0e@ornl.gov

Steps to run the systems:

Step1: Use geninit_config/initialize_dirs_for_run.py to generate initial configuration for the lignin melts at high temperature. This requires all the files from SPRInG in the same directory and supp_initdirs.py
Step2: Use myruns/myruns/tg_input.py and set the required temperature range to generate lignin melt simulations at different temperature.
Step3: Use gmx_functions.py to compute Rg/MSD/segmentalMSD
Step4: Use plotdata.py to plot data and to compute Rg-N scaling and Tg values.

