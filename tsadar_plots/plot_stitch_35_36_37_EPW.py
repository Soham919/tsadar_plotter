import numpy as np
import matplotlib.pyplot as plt
import cv2
import os 
from pathlib import Path
import re
import pandas as pd
from TS_aux import load_file
from TS_aux import shot_data


def plot_stitch_35_36_37_EPW(data):
    ### PLOT STITCHED H2 SHOCKS ###
    im1 = data['92536']['92536_9_nersc.csv']
    im2 = data['92537']['92537_8_nersc.csv']
    #im3 = data['92537']['92537_5_nersc.csv']
    magI = 2.87
    r1 = shot_data.loc[92536,"pointing"] - im1["Radius (\\mum)"]/10**3  # convert to mm and shift 
    r2 = shot_data.loc[92537,"pointing"] - im2["Radius (\\mum)"]/10**3 
    #r3 = shot_data.loc[92537,"pointing"] - im3["Radius (\\mum)"]/10**3
    ## Ti, Va, fract plot for all 3 species 
    plt.rcParams.update({
    #    "font.family": "serif",             # pick font family
    #    "font.serif": ["Arial"],  # pick specific font
        "font.size": 17,                     # set size
        "axes.labelsize": 17,    # font size for axis labels
        "axes.titlesize": 17,     # font size for titles
        # Tick direction
        "xtick.direction": "in",
        "ytick.direction": "in",
        # Tick width
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        # Tick length
        "xtick.major.size": 3.5,
        "ytick.major.size": 3.5,
    })

    fig1, ax1  = plt.subplots(figsize=(16,13))
    # Error parameters

    window = 3  # rolling window
    f1_1 = im1["Te_electron"].rolling(window=window).mean()
    f2_1 = im1["ne_electron"].rolling(window=window).mean()
    f1_2 = im2["Te_electron"].rolling(window=window).mean()
    f2_2 = im2["ne_electron"].rolling(window=window).mean()

    std_f1_1 = im1["Te_electron"].rolling(window=window).std()
    std_f2_1 = im1["ne_electron"].rolling(window=window).std()
    std_f1_2 = im2["Te_electron"].rolling(window=window).std()
    std_f2_2 = im2["ne_electron"].rolling(window=window).std()


    fig1, ax1  = plt.subplots(figsize=(7.5,5.9))
    c1 = "tab:red"
    c2 = "tab:blue"
    l1, = ax1.plot(r1,f1_1,color=c1,linewidth=2,linestyle='-',label=r'$T_{e}$')
    l2, = ax1.plot(r2,f1_2,color=c1,linewidth=2,linestyle='--',label=r'$T_{e}$')
    ax1.fill_between(r1,f1_1-std_f1_1,f1_1+std_f1_1,
                    color='LightCoral',        # fill color
                    alpha=0.2,           # transparency
                    edgecolor='IndianRed',    # outline color
                    linewidth=1.5,       # outline width
                    linestyle='-',      # outline style (optional)
                    label='±1σ band')
    
    ax1.fill_between(r2,f1_2-std_f1_2,f1_2+std_f1_2,
                    color='LightCoral',        # fill color
                    alpha=0.2,           # transparency
                    edgecolor='IndianRed',    # outline color
                    linewidth=1.5,       # outline width
                    linestyle='-',      # outline style (optional)
                    label='±1σ band')

    #ax1[2].plot(r, lp["Va_ion-2"],color='mediumblue',linewidth=2,linestyle='-')
    #ax1[2].plot(r, bv,color='k',linewidth=2,linestyle='--')
    ax1.set_title('Electron density and temperature')
    ax1.set_xlabel(r"$x (mm)$")
    #ax1.set_xlim([3.55, 4.75])
    #ax1[1,1].set_ylim([0.22, 0.4])
    ax1.invert_xaxis()
    ax1.set_ylabel(r"$T_{e} (keV)$", color = c1)
    ax1.tick_params(axis='y', colors=c1)
    ax1.spines['left'].set_color(c1)
    ax1.spines['left'].set_linewidth(1.5)
    ax1.spines['right'].set_visible(False)
    ax1.grid(True)
    #ax1[1,1].legend()

    ax2 = ax1.twinx()
    l3, = ax2.plot(r1,f2_1,color=c2,linewidth=2,linestyle='-',label=r'$n_{e}$')
    l4, = ax2.plot(r2,f2_2,color=c2,linewidth=2,linestyle='--',label=r'$n_{e}$')
    ax2.fill_between(r1,f2_1-std_f2_1,f2_1+std_f2_1,color='DodgerBlue',        # fill color
                    alpha=0.2,           # transparency
                    edgecolor='RoyalBlue',    # outline color
                    linewidth=1.5,       # outline width
                    linestyle='-',      # outline style (optional)
                    label='±1σ band')
    ax2.fill_between(r2,f2_2-std_f2_2,f2_2+std_f2_2,color='DodgerBlue',        # fill color
                    alpha=0.2,           # transparency
                    edgecolor='RoyalBlue',    # outline color
                    linewidth=1.5,       # outline width
                    linestyle='-',      # outline style (optional)
                    label='±1σ band')

    ax2.set_ylabel(r"$n_{e} (\times 10^{20} cm^{-3})$", color = c2)
    ax2.tick_params(axis='y', colors=c2)
    ax2.spines['right'].set_color(c2)
    ax2.spines['right'].set_linewidth(1.5)
    ax2.spines['left'].set_visible(False)
    #ax2.legend()

    lines = [l1, l3]
    labels = [line.get_label() for line in lines]
    ax1.legend(lines, labels, loc='upper right', frameon=False)








   
    