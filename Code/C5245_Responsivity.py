# -*- coding: utf-8 -*-
"""
Created on Wed Dec 13 14:36:56 2023

@author: pwilson3
"""

#%% PACKAGES AND GLOBAL SETTINGS ################################################

from cycler import cycler
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy import interpolate
import numpy as np
import csv
from datetime import datetime
from matplotlib.lines import Line2D
from scipy import constants as const
import os
import sys

# Add the paths to where you keep the measurements (.sciv files) and to where you
# are keeping the PPC_data_collector.py and plot_PPC_IV_functions.py files.
sys.path.append('C:\\Users\\louis\\OneDrive - University of Ottawa\\uOttawa\\2023-24\\Hiver2024\\Sunab pt 2\\Data\\') #!!!
sys.path.append('C:\\Users\\louis\\Sunlab Analysis Local\\Two-junction\\Getting_started_files') #!!!
from Sample_Names_and_Locations import *
import PPC_data_collector as PPC
import plot_PPC_IV_functions as plt_PPC

#%% DATA LOCATIONS ##############################################################

# Path to where the data is located. Same as in the import section.
directory = 'C:\\Users\\louis\\OneDrive - University of Ottawa\\uOttawa\\2023-24\\Hiver2024\\Sunab pt 2\\Data\\' #!!!

# Laser power measurements
# Path to where the 1520 nm laser power file is being kept.
power_folder = 'C:\\Users\\louis\\OneDrive - University of Ottawa\\uOttawa\\2023-24\\Hiver2024\\Sunab pt 2\\Data\\' #!!!
power_file = '2023-10-02-Compiled_Power_Measurements_for_LTT_1520nm.xlsx'
power_table = pd.read_excel(power_folder+power_file, 
                            sheet_name = 'Values_for_Current_Use',
                            header=2, usecols='B, C, D, E, F, G, H, I, J')

power_table = power_table.set_index('Current target (A)')
power_table = power_table.dropna()


# This file has estimates for the 1310 nm laser power.
#-------------------------------
# power_file2 = 'Power_Measurements_for_1310nm.xlsx'
# power_table_1310nm = pd.read_excel(power_folder+power_file2, 
#                                             header=2, usecols='B, D',
#                                             sheet_name = 'Sheet1')
# power_table_1310nm = power_table_1310nm.set_index('Unnamed: 1')
# power_table_1310nm = power_table_1310nm.dropna()
#-------------------------------
# This file has the thermopile power measurements
# This is probably not needed for most instances but can be good for comparing 
# calibrated vs. thermopile power measurements.
power_table_orig = pd.read_excel(power_folder+power_file, 
                            sheet_name = 'Calibrations',
                            header=2, usecols='B, C, D, E, F, G, H, I, J')
power_table_orig = power_table_orig.set_index('Current target (A)')
power_table_orig = power_table_orig.dropna()

#%% GET DATA, ORGANIZE, AND CLEANUP ############################################

# Samples you want to look at.
# These lists of folders can be defined in the Sample_Names_and_Locations file
# or elsewhere. This just puts them all together.
folderlist = folder_C5245 + folder_C5246 + folder_C5247 + folder_C5195

# Put all the data into an initial dictionary. Choose what type of dictionary
# you want (DarkIV, LightIV, FixedCurrent, FixedVoltage, LightBias) or combine
# dictionaries by adding multiple file types to one list.

DarkIV_data = PPC.PPC_data_collector(['DarkIV'], directory, folderlist)

LightIV_data = PPC.PPC_data_collector(['LightIV'], directory, folderlist)

FixedVoltage_data = PPC.PPC_data_collector(['FixedVoltage'],directory,folderlist)
#Get the fixed voltage data
FixedVoltage_data.get_powers_updated(power_table,
                                alternate_table = power_table_orig,
                                filters = 'ND1 + ND2', 
                                keyND1ND2 = 'ND1+ND2', 
                                keyUncertainty = 'Uncertainty',
                                AlternatekeyUncertainty = 'Uncertainty')
FixedVoltage_data.get_powers_updated(power_table, filters = 'ND2', 
                                alternate_table = power_table_orig,
                                keyND2 = 'ND2', 
                                keyUncertainty = 'Uncertainty',
                                AlternatekeyUncertainty = 'Uncertainty')
FixedVoltage_data.get_powers_updated(power_table, filters = 'ND1', 
                                alternate_table = power_table_orig,
                                keyND1 = 'ND1', 
                                keyUncertainty = 'Uncertainty',
                                AlternatekeyUncertainty = 'Uncertainty')
FixedVoltage_data.get_powers_updated(power_table, filters = 'No Filter', 
                                alternate_table = power_table_orig,
                                keyNoFilter = 'Power No Filter (W)', 
                                keyUncertainty = 'Uncertainty',
                                AlternatekeyUncertainty = 'Uncertainty')

# Add the 1550 nm laser powers to the dictionary.
LightIV_data.get_powers_updated(power_table,
                                alternate_table = power_table_orig,
                                filters = 'ND1 + ND2', 
                                keyND1ND2 = 'ND1+ND2', 
                                keyUncertainty = 'Uncertainty',
                                AlternatekeyUncertainty = 'Uncertainty')
LightIV_data.get_powers_updated(power_table, filters = 'ND2', 
                                alternate_table = power_table_orig,
                                keyND2 = 'ND2', 
                                keyUncertainty = 'Uncertainty',
                                AlternatekeyUncertainty = 'Uncertainty')
LightIV_data.get_powers_updated(power_table, filters = 'ND1', 
                                alternate_table = power_table_orig,
                                keyND1 = 'ND1', 
                                keyUncertainty = 'Uncertainty',
                                AlternatekeyUncertainty = 'Uncertainty')
LightIV_data.get_powers_updated(power_table, filters = 'No Filter', 
                                alternate_table = power_table_orig,
                                keyNoFilter = 'Power No Filter (W)', 
                                keyUncertainty = 'Uncertainty',
                                AlternatekeyUncertainty = 'Uncertainty')

# Add the 1310 nm power to the dictioanry.
# LightIV_data.get_powers_updated(power_table_1310nm, LaserWL = '1310nm',
#                                 keyNoFilter1310= 'Power No Filter (W)')

# Add the 1550 nm + 1310 nm bias light conditions to the dictionary.
# LightBias_data.get_powers_updated(power_table, 
#                         alternate_table = power_table_1310nm,
#                         filters = 'No Filter', LaserWL = '1310nm+1550nm',
#                         keyNoFilter = 'Power No Filter (W)',
#                         keyNoFilter1310= 'Power No Filter (W)')
# LightBias_data.get_powers_updated(power_table, 
#                         alternate_table = power_table_1310nm,
#                         filters = 'ND1', LaserWL = '1310nm+1550nm',
#                         keyND1 = 'ND1',
#                         keyNoFilter1310= 'Power No Filter (W)')
# LightBias_data.get_powers_updated(power_table, 
#                         alternate_table = power_table_1310nm,
#                         filters = 'ND2', LaserWL = '1310nm+1550nm',
#                         keyND2 = 'ND2',
#                         keyNoFilter1310= 'Power No Filter (W)')
# LightBias_data.get_powers_updated(power_table, 
#                         alternate_table = power_table_1310nm,
#                         filters = 'ND1 + ND2', LaserWL = '1310nm+1550nm',
#                         keyND1ND2 = 'ND1+ND2',
#                         keyNoFilter1310= 'Power No Filter (W)')

# Calculate some initial values (efficiency, etc.) and add them to the 
# dictionary. Option to use a different power value (for example to compare 
# with calibrated vs. uncalibrated values.)  
LightIV_data.add_to_dict('initial', power_label = 'Power (W)', 
                         alternate_power_label = 'Power (LTT tuned) (W)')
                         
# Some other optional things you can calculate but may not need to.
LightIV_data.add_to_dict('fit_Voc_slope', Atot = 0.054)
LightIV_data.add_to_dict('fit_Jsc_slope', Atot = 0.054)
#LightBias_data.add_to_dict('initial', power_label = 'Power (W)')

Atot = 0.054 # (cm^2)
Atot2 = 1 # (cm^2)

voltageDic = FixedVoltage_data.dictionary
lightDic= LightIV_data.dictionary

plt.close('all')

#%% PLOT GLOBAL SETTINGS ########################################################
color = sns.color_palette('Set1')
colours = (cycler(color=sns.color_palette('Set1')))
color_grad = sns.color_palette('rocket',  n_colors = 130)
color_grad2 = sns.color_palette('Blues_r', n_colors = 130)
# colors = (cycler(color=sns.color_palette('rocket',  n_colors = 48)))
#
sample = ['C5245-X7Y0'] # 'C5246-X12Y1', 'C5247-X7Y6'
sample_labels = ['C5245-X7Y0'] # 'C5246-X12Y1', 'C5247-X7Y6' #For plotting
dates = [None, None, None, None, None, None, None, None, None, None, None]
colour_C5195 = color[0]
colour_C5245 = color[1]
colour_C5246 = color[2]
colour_C5247 = color[3]
colours = [colour_C5245,colour_C5245, colour_C5246, colour_C5246, colour_C5247, colour_C5247, colour_C5247, colour_C5247, colour_C5195, colour_C5195]
markers =['o','o','o','o','o','o','o','o','o','o','o','o','o','o','o']#['.', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X'] #['None','None','None','None','None','None','None','None','None','None','None','None','None','None','None']#['.', 'o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'P', 'X']
linestyles = ['solid', 'dashed','solid','dashed', 'solid','dotted','dashed','dashdot','solid','dashed']
linewidth = 2
markersize = 3
legend_fontsize = 15

sns.set_theme(context='poster', style="ticks",
              rc={'lines.markeredgecolor': 'k', 'axes.edgecolor': 'k',
                  'xtick.direction': 'in', 'ytick.direction': 'in',
                  'xtick.top': True, 'ytick.right': True,
                  'lines.markersize': markersize})

#%% PLOT SPECTRAL RESPONSE VS IRRADIANCE####################################################

fig1, ax1 = plt.subplots(1, 1)
fig1.set_size_inches(11, 9)
# ax1.set_prop_cycle(colors)
# ax1.axhline(0, c = 'k')
# ax1.axvline(0, c = 'k')

#sample = ['C5245-X7Y0','C5245-X3Y1', 'C5245-X6Y1', 'C5246-X12Y1','C5246-X8Y1', 'C5247-X7Y6', 'C5247-X7Y5', 'C5247-X5Y5', 'C5247-X5Y4'] # 'C5246-X12Y1', 'C5247-X7Y6'
#sample_labels = ['C5245-X7Y0','C5245-X3Y1', 'C5245-X6Y1', 'C5246-X12Y1','C5246-X8Y1', 'C5247-X7Y6', 'C5247-X7Y5', 'C5247-X5Y5', 'C5247-X5Y4'] # 'C5246-X12Y1', 'C5247-X7Y6' #For plotting
#dates = [None, None, None, None, None, None, None, None, None]
#colours = [color[0], color[1], color[2], color[3], color[4], color[5], color[6], color[7], color[8]]
#markers = ['o', 's', '*', '^', 'P', 'o', 's', '*', '^']

for i in range(len(sample)):
    # This grab_data function goes into the data and grabs certain parameters 
    # you specifiy. In this case the power and efficiency.
    # Select the dictionary you want to pull from, the  focus lens position, 
    # the filter(s) you want, and the laser current(s) you want to pull
    x1, y1 = plt_PPC.grab_data(FixedVoltage_data.dictionary, sample[i], ['NaN'], ['ND1 + ND2','ND2','ND1','No Filter'],['6.0','8.0','10.0','12.0','14.0','16.0','18.0', '6.5', '7.0'],
                          'Power (W)', 'Current (A)',
                          # pull out data only from a specific date (optional)
                          date =  dates[i],
                          # pull out data only with a specific wavelength (1550nm or 1310nm)
                          LaserWL = '1550nm',
                          # pull out data only with a specific voltage sweep direction
                          direction = 'Forward',
                          # x and y factor scales (in this case turn power into irradiance)
                          xfactor=1/Atot, yfactor=-1, 
                          # the below won't do anything here
                          m='^', col=color[3],
                          lab='nolabel', l='')
    n = 0
    y_prime = [None]*len(y1)
    for y in y1:
        temp = y[-1]
        y_prime[n] = temp
        n += 1
    #resp = np.divide(y2,x1)
    ax1.plot(x1, np.divide(y_prime,x1), c = colours[i], marker = markers[i], ls = linestyles[i], lw = linewidth,
             label = "25$\degree$C")
    x2, y2 = plt_PPC.grab_data(FixedVoltage_data.dictionary, sample[i], ['NaN'], ['ND1 + ND2','ND2','ND1','No Filter'],['6.0 (1)','8.0 (1)','10.0 (1)','12.0 (1)','14.0 (1)','16.0 (1)','18.0 (1)', '6.5 (1)', '7.0 (1)'],
                          'Power (W)', 'Current (A)',
                          # pull out data only from a specific date (optional)
                          date =  dates[i],
                          # pull out data only with a specific wavelength (1550nm or 1310nm)
                          LaserWL = '1550nm',
                          # pull out data only with a specific voltage sweep direction
                          direction = 'Forward',
                          # x and y factor scales (in this case turn power into irradiance)
                          xfactor=1/Atot, yfactor=-1, 
                          # the below won't do anything here
                          m='^', col=color[3],
                          lab='nolabel', l='')
    n = 0
    y_prime = [None]*len(y2)
    for y in y2:
        temp = y[-1]
        y_prime[n] = temp
        n += 1
    #resp = np.divide(y2,x1)
    ax1.plot(x2, np.divide(y_prime,x2), c = colours[i+4], marker = markers[i+4], ls = linestyles[i], lw = linewidth,
             label = "23$\degree$C")
    x3, y3 = plt_PPC.grab_data(FixedVoltage_data.dictionary, sample[i], ['NaN'], ['ND1 + ND2','ND2','ND1','No Filter'],['6.0 (2)','8.0 (2)','10.0 (2)','12.0 (2)','14.0 (2)','16.0 (2)','18.0 (2)', '6.5 (2)', '7.0 (2)'],
                          'Power (W)', 'Current (A)',
                          # pull out data only from a specific date (optional)
                          date =  dates[i],
                          # pull out data only with a specific wavelength (1550nm or 1310nm)
                          LaserWL = '1550nm',
                          # pull out data only with a specific voltage sweep direction
                          direction = 'Forward',
                          # x and y factor scales (in this case turn power into irradiance)
                          xfactor=1/Atot, yfactor=-1, 
                          # the below won't do anything here
                          m='^', col=color[3],
                          lab='nolabel', l='')
    n = 0
    y_prime = [None]*len(y2)
    for y in y3:
        temp = y[-1]
        y_prime[n] = temp
        n += 1
    #resp = np.divide(y2,x1)
    ax1.plot(x3, np.divide(y_prime,x3), c = 'orange', marker = markers[i+7], ls = linestyles[i], lw = linewidth,
             label = "21$\degree$C")
    x, y = plt_PPC.grab_data(LightIV_data.dictionary, sample[i], ['18mm'], ['all'], ['all'],
                          'Power (W)', 'Isc (A)',
                          # pull out data only from a specific date (optional)
                          date =  dates[i],
                          # pull out data only with a specific wavelength (1550nm or 1310nm)
                          LaserWL = '1550nm',
                          # pull out data only with a specific voltage sweep direction
                          direction = 'Forward',
                          # x and y factor scales (in this case turn power into irradiance)
                          xfactor=1/Atot, yfactor=-1, 
                          # the below won't do anything here
                          m='^', col=color[3],
                          lab='nolabel', l='')
    y = abs(np.array(y))
    ax1.plot(x, np.divide(y,x), c = colours[i+3], marker = markers[i+3], ls = linestyles[i], lw = linewidth,
             label = sample_labels[i])
    

ax1.set_xscale('log')
ax1.grid(which='both')
ax1.set_xlabel('Irradiance (W/cm$^{2}$)')
ax1.set_ylabel('Responsivity (A/W)')
ax1.set_title("Responsivity of various 2J PPCs")
fig1.legend(framealpha=1, loc = "upper left", fontsize = legend_fontsize).set_draggable(True)
plt.show()