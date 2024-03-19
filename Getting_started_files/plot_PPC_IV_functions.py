8.# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 16:07:32 2023

@author: pwilson3
"""

import pandas as pd
from os import listdir, path
import matplotlib.pyplot as plt
import seaborn as sns
from cycler import cycler
import seaborn as sns
import numpy as np
from scipy.optimize import curve_fit

# FUNCTIONS

def grab_data(dataset, name, positions, filters, currents, xparam, yparam,
                  direction = None, date = None, laser_wavl = None, 
                  xfactor = 1, yfactor = 1, 
                  l = '-', m = None, col = None, 
                  lab =  None, LaserWL = '1550nm', script_name = 'LightIV',
                  LaserB = None):
    """
    This grabs data in the exat same way the plot_params does except it does't 
    plot the data it just returns the x and y values'

    Parameters
    ----------
    dataset : DICTIONARY THAT CONTAINS ALL THE DATA (ALL_DATA, ETC.).
    name : SAMPLE NAME (EX. 'C5031-X0Y5')
    positions : FILTER POSITION AS STRING ('18mm')
    filters : FILTERS IN PLACE
    currents : LASER TARGET CURRENT
    xparam : X-AXIS DATA 
    yparam : Y-AXIS DATA
    direction : TYPE, optional
        DESCRIPTION. The default is None.
    date : TYPE, optional
        DESCRIPTION. The default is None.
    xfactor : TYPE, optional
        DESCRIPTION. The default is 1.
    yfactor : TYPE, optional
        DESCRIPTION. The default is 1.
    l : TYPE, optional
        DESCRIPTION. The default is '-'.
    m : TYPE, optional
        DESCRIPTION. The default is None.
    col : TYPE, optional
        DESCRIPTION. The default is None.
    lab : TYPE, optional
        DESCRIPTION. The default is None.

    Returns
    -------
    None.

    """
    x = []
    y = []
    pos = positions
    if pos == ['all']:
        pos = dataset[name].keys()
    for p in pos:
        filt = filters
        if filt == ['all']:
            filt = dataset[name][p].keys()
        for f in filt:
            curr = currents
            if direction != None or date != None or laser_wavl != None:
                if curr == ['all']:
                    curr = dataset[name][p][f].keys()
                #if curr[0][0] 

                if script_name == 'LightIV':
                    print("Fails at", name)
                    curr = [x for x in curr if 
                        (dataset[name][p][f][x]['Direction'] == direction or direction == None)
                        and (dataset[name][p][f][x]['Date'] == date or date == None)
                        and (dataset[name][p][f][x]['Laser Wavelength'] == LaserWL)]
                elif script_name == 'LightBias':
                    curr = [x for x in curr if 
                            (dataset[name][p][f][x]['Direction'] == direction or direction == None)
                            and (dataset[name][p][f][x]['Date'] == date or date == None)
                            and (dataset[name][p][f][x]['Laser B Current (A)'] == LaserB)]
                    # print('!!!:', [(dataset[name][p][f][x]['Direction'] == direction or direction == None)
                            # and (dataset[name][p][f][x]['Date'] == date or date == None)
                            # and (dataset[name][p][f][x]['Laser B Current (A)'] == LaserB) for x in curr])
                else:
                    try:
                        curr = [x for x in curr if 
                                (dataset[name][p][f][x]['Date'] == date or date == None)
                                and (dataset[name][p][f][x]['Laser Wavelength'] == LaserWL)]
                    except:
                        print("Could not find a current...")
                try:
                    curr = [x for x in curr if 
                            (dataset[name][p][f][x]['Direction'] == direction or direction == None)
                            and (dataset[name][p][f][x]['Date'] == date or date == None)
                            and (dataset[name][p][f][x]['Laser Wavelength'] == laser_wavl or laser_wavl == None)]
                except:
                    print("Whoops")
            else:
                if curr == ['all']:
                    try:
                        curr = dataset[name][p][f].keys()
                    except KeyError as err:
                        print('key error in plot_params:', err)
                        print(name, p, f, curr)
                        

            for c in curr:
                if xparam != 'Laser Current':
                    try:
                        y.append(yfactor*dataset[name][p][f][c][yparam])
                        x.append(xfactor*dataset[name][p][f][c][xparam])
                        # print(name, p, f, c, xparam, yparam)
                    except KeyError as err:
                        print('key error in plot_params:', err)
                        print(name, p, f, c)
                else:
                    try:
                        y.append(yfactor*dataset[name][p][f][c][yparam])
                        x.append(xfactor*float(c.split('(')[0]))
                        # print(x)
                        # print(y)
                        # print(name, p, f, c, xparam, yparam)
                    except KeyError as err:
                        print('key error in plot_params:', err)
                        print(name, p, f, c)
    return x, y


def plot_PPC_data(ax, dataset, name, positions, filters, currents, xparam, yparam,
                  direction = None, date = None, laser_wavl = None, xfactor = 1, yfactor = 1, l = '-', m = None, col = None, 
                  lab =  None, LaserWL = '1550nm', script_name = 'LightIV',
                  LaserB = None, slice_index = [None, None], shift_zero = False,
                  slice_type = None):
    pos = positions
    if pos == ['all']:
        pos = dataset[name].keys()
    for p in pos:
        filt = filters
        if filt == ['all']:
            filt = dataset[name][p].keys()
        for f in filt:
            curr = currents
            if direction != None or date != None or laser_wavl != None:
                curr = currents
                if curr == ['all']:
                    curr = dataset[name][p][f].keys()
                if script_name == 'LightIV':
                    curr = [x for x in curr if 
                            (dataset[name][p][f][x]['Direction'] == direction or direction == None)
                            and (dataset[name][p][f][x]['Date'] == date or date == None)
                            and (dataset[name][p][f][x]['Laser Wavelength'] == LaserWL)]
                elif script_name == 'LightBias':
                    curr = [x for x in curr if 
                            (dataset[name][p][f][x]['Direction'] == direction or direction == None)
                            and (dataset[name][p][f][x]['Date'] == date or date == None)
                            and (dataset[name][p][f][x]['Laser Wavelength'] == LaserWL)
                            and (dataset[name][p][f][x]['Laser B Current (A)'] == LaserB)]
                else:
                    curr = [x for x in curr if 
                            (dataset[name][p][f][x]['Date'] == date or date == None)
                            and (dataset[name][p][f][x]['Laser Wavelength'] == LaserWL)]
                # print(curr)
                if type(col) == sns.palettes._ColorPalette:
                    max_curr = float(curr[-1].split('(')[0])
                    # print(max_curr)
                    # print(len(col))
                    col_pal = [col[int(100*(float(c.split('(')[0]))/max_curr)] for c in curr]
                    # print([int(100*(float(c.split('(')[0]))/max_curr) for c in curr])
                    ax.set_prop_cycle(cycler("color", sns.color_palette(col_pal)))
                    col = None
                for c in curr:
                    if slice_type == 'from_beam_block_file':
                        slice_index = [dataset[name][p][f][c]['Unblocked Index'], None]
                        # print(slice_index)
                    shift = 0
                    if shift_zero == True:
                        shift = dataset[name][p][f][c]['Time Blocked (s)']
                    try:
                        ax.plot(xfactor*dataset[name][p][f][c][xparam][slice_index[0]:slice_index[1]] - shift, 
                                 yfactor*dataset[name][p][f][c][yparam][slice_index[0]:slice_index[1]],
                                 ls = l, marker = m, color = col, label = lab)
                        
                    except KeyError as err:
                        print('key error:', err)
                        print(name, p, f, c)
            else:
                curr = currents
                if curr == ['all']:
                    curr = dataset[name][p][f].keys()
                # print(curr)
                if type(col) == sns.palettes._ColorPalette:
                    max_curr = float(curr[-1].split('(')[0])
                    # print(max_curr)
                    # print(len(col))
                    col_pal = [col[int(100*(float(c.split('(')[0]))/max_curr)] for c in curr]
                    # print([int(100*(float(c.split('(')[0]))/max_curr) for c in curr])
                    ax.set_prop_cycle(cycler("color", sns.color_palette(col_pal)))
                    col = None
                for c in curr:
                    if slice_type == 'from_beam_block_file':
                        slice_index = [dataset[name][p][f][c]['Unblocked Index'], None]
                        # print(slice_index)
                    shift = 0
                    if shift_zero == True:
                        shift = dataset[name][p][f][c]['Time Blocked (s)']
                    try:
                        ax.plot(xfactor*np.array(dataset[name][p][f][c][xparam][slice_index[0]:slice_index[1]] - shift), 
                                 yfactor*np.array(dataset[name][p][f][c][yparam][slice_index[0]:slice_index[1]]),
                                 ls = l, marker = m, color = col, label = lab)
                    except KeyError as err:
                        print('key error:', err)
                        print(name, p, f, c)

def plot_params(ax, dataset, name, positions, filters, currents, xparam, yparam,
                  direction = None, date = None, xfactor = 1, yfactor = 1, l = '-', m = None, col = None, 
                  lab =  None, a  = None, LaserWL = '1550nm', script_name = 'LightIV',
                  LaserB= None, x_error = None, y_error = None):
    x = []
    y = []
    x_err = []
    y_err = []
    pos = positions
    if pos == ['all']:
        pos = dataset[name].keys()
    # print(positions)
    for p in pos:
        filt = filters
        if filt == ['all']:
            filt = dataset[name][p].keys()
        # print(filters)
        for f in filt:
            curr = currents
            if (direction != None or date != None or LaserB != None):
                if curr == ['all']:
                    curr = dataset[name][p][f].keys()
                if script_name == 'LightIV':
                    curr = [x for x in curr if 
                            (dataset[name][p][f][x]['Direction'] == direction or direction == None)
                            and (dataset[name][p][f][x]['Date'] == date or date == None)
                            and (dataset[name][p][f][x]['Laser Wavelength'] == LaserWL)]
                elif script_name == 'LightBias':
                    curr = [x for x in curr if 
                            (dataset[name][p][f][x]['Direction'] == direction or direction == None)
                            and (dataset[name][p][f][x]['Date'] == date or date == None)
                            and (dataset[name][p][f][x]['Laser Wavelength'] == LaserWL)
                            and (dataset[name][p][f][x]['Laser B Current (A)'] == LaserB)]
                else:
                    curr = [x for x in curr if 
                            (dataset[name][p][f][x]['Date'] == date or date == None)
                            and (dataset[name][p][f][x]['Laser Wavelength'] == LaserWL)]
            else:
                if curr == ['all']:
                    curr = dataset[name][p][f].keys()
            # print(currents)
            for c in curr:
                if xparam != 'Laser Current':
                    try:
                        # print(x)
                        # print(y)
                        y.append(yfactor*dataset[name][p][f][c][yparam])
                        x.append(xfactor*dataset[name][p][f][c][xparam])
                        
                        if x_error != None:
                            x_err.append(xfactor*dataset[name][p][f][c][x_error])
                        if y_error != None:
                            y_err.append(yfactor*dataset[name][p][f][c][y_error])
                        # print(name, p, f, c, xparam, yparam)
                    except KeyError as err:
                        print('key error in plot_params:', err)
                        print(name, p, f, c)
                else:
                    try:
                        y.append(yfactor*dataset[name][p][f][c][yparam])
                        x.append(xfactor*float(c.split('(')[0]))
                        print(x)
                        print(y)
                        # print(name, p, f, c, xparam, yparam)
                    except KeyError as err:
                        print('key error in plot_params:', err)
                        print(name, p, f, c)
    ax.plot(x, y,
             ls = l, marker = m, color = col, label = lab, alpha = a)
    if x_error != None or y_error != None:
        if len(x_err) == 0:
            x_err = None
        if len(y_err) == 0:
            y_err = None
        lines, caps, bars = ax.errorbar(x, y, xerr = x_err, yerr = y_err, ls = '', ecolor = 'k', 
                    capsize = 3.0, capthick = 2.0, elinewidth = 2)

def linear_fit(x, m, b):
    return m*x + b

def linear_fit_yint_zero(m, x):
    return m[0]*x

def poly(x, a, b, c):
    return a + b*x + c*x*x

def fit_Pmax(dataset, name, positions, filters, currents, xparam, yparam):
    """
    Finds Pmax with a fit

    Parameters
    ----------
    start : P_MAX TAKE FROM IV FILE PARAMETERS

    Returns
    -------
    None.

    """
    
    x, y = grab_data(dataset, name, positions, filters, currents, xparam, yparam)
    x = x[0] # this is necessary because grab_data was based off the plot_params functions which puts things in a list
    y = y[0]
    # print(x)
    # print(y)
    max_index = np.where(y == min(y))[0][0]
    # print(max_index)
    slicex = x[max_index - 1:max_index + 2] # grab a couple of data points around the maximum
    slicey = y[max_index - 1:max_index + 2]
    # print(slicex)
    # print(slicey)
    params = curve_fit(poly, slicex, slicey)
    newx = np.linspace(slicex[0], slicex[-1], 1000)
    newy = poly(newx, params[0][0], params[0][1], params[0][2])
    new_max = min(newy)
    new_maxx = newx[np.where(newy == min(newy))][0]
    return new_max, new_maxx, newx, newy

def fit_Isc(dataset, name, positions, filters, currents, xparam, yparam):
    """
    Finds Pmax with a fit

    Parameters
    ----------
    start : P_MAX TAKE FROM IV FILE PARAMETERS

    Returns
    -------
    None.

    """
    
    x, y = grab_data(dataset, name, positions, filters, currents, xparam, yparam)
    x = x[0] # this is necessary because grab_data was based off the plot_params functions which puts things in a list
    y = y[0]
    # print(x)
    # print(y)
    zero_index = np.where(np.absolute(x) == min(np.absolute(x)))[0][0] # find point closest to V = 0
    # print(zero_index)
    slicex = x[zero_index - 1:zero_index + 2] # grab a couple of data points around the maximum
    slicey = y[zero_index - 1:zero_index + 2]
    # print(slicex)
    # print(slicey)
    params = curve_fit(linear_fit, slicex, slicey)
    newx = np.linspace(slicex[0], slicex[-1], 1000)
    newy = linear_fit(newx, params[0][0], params[0][1])
    new_zero = newx[np.where(np.absolute(newx) == min(np.absolute(newx)))][0]
    new_Isc = newy[np.where(np.absolute(newx) == min(np.absolute(newx)))][0]
    return new_Isc, new_zero, newx, newy

def fit_Voc_slope(dataset, name, positions, filters, currents, xparam, yparam,
                  yfactor):
    """
    Finds the slope of the current at Voc (approximate)

    Parameters
    ----------
    start : dictionary parameters for the sample

    Returns
    -------
    None.

    """
    
    x, y = grab_data(dataset, name, positions, filters, currents, xparam, yparam,
                     yfactor = -1*yfactor) # xparam, yparam should be voltage and current (current x -1 so the slope will be negative)
    x = x[0] # this is necessary because grab_data was based off the plot_params functions which puts things in a list
    y = y[0]
    # print(x)
    # print(y)
    max_index = np.where(np.absolute(y) == min(np.absolute(y)))[0][0] # find the voltage index for value close to Voc
    # print(max_index)
    if y[max_index] > 0:
        slicex = x[max_index:max_index + 2] # grab a couple of data points around Voc
        slicey = y[max_index:max_index + 2]
    else:
        slicex = x[max_index - 1:max_index+1] # grab a couple of data points around Voc
        slicey = y[max_index - 1:max_index+1]
    # print(slicex)
    print(slicey)
    print(name, positions, filters, currents)
    params = curve_fit(linear_fit, slicex, slicey)
    slope_at_Voc = params[0][0]
    fit_intercept = params[0][1]
    return slope_at_Voc, fit_intercept

def fit_Jsc_slope(dataset, name, positions, filters, currents, xparam, yparam,
                  yfactor):
    """
    Finds the slope of the current at Jsc (approximate)

    Parameters
    ----------
    start : dictionary parameters for the sample

    Returns
    -------
    None.

    """
    
    x, y = grab_data(dataset, name, positions, filters, currents, xparam, yparam,
                     yfactor = -1*yfactor) # xparam, yparam should be voltage and current (current x -1 so the slope will be negative)
    x = x[0] # this is necessary because grab_data was based off the plot_params functions which puts things in a list
    y = y[0]
    # print(x)
    # print(y)
    zero_index = np.where(np.absolute(x) == min(np.absolute(x)))[0][0] # find the voltage index for value close to Jsc
    # print(max_index)
    slicex = x[zero_index-1:zero_index + 4] # grab a couple of data points around Voc
    slicey = y[zero_index-1:zero_index + 4]

    # print(slicex)
    print(slicey)
    print(name, positions, filters, currents)
    params = curve_fit(linear_fit, slicex, slicey)
    slope_at_Jsc = params[0][0]
    fit_intercept = params[0][1]
    return slope_at_Jsc, fit_intercept