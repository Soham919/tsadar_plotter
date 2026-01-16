import numpy as np
import matplotlib.pyplot as plt
import cv2
import os 
from pathlib import Path
import re
import pandas as pd
from TS_aux import load_file
from TS_aux import shot_data
from Plasma_func import f_rho

def plot_92537(data):
    ### PLOT LEARNED PARAMETERS for H2 AND HE SHOCKS ###
    im1 = data['92537']['92537_5_nersc.csv']
    epw = data['92537']['92537_8_nersc.csv']

    # Convert ion fractions to number denisties




    r = shot_data.loc[92537,"pointing"] - im1["Radius (\\mum)"]/10**3
    bv = (r*10**2)/shot_data.loc[92537,"timing"]
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

    fig1, ax1  = plt.subplots(2,2,figsize=(16,13))

    window = 2  # rolling window
    f1_1 = im1["fract_ion-1"].rolling(window=window).mean()
    f2_1 = im1["fract_ion-2"].rolling(window=window).mean()
    #f3_1 = im1["fract_ion-3"].rolling(window=window).mean()
    std_f1_1 = im1["fract_ion-1"].rolling(window=window).std()
    std_f2_1 = im1["fract_ion-2"].rolling(window=window).std()
    #std_f3_1 = im1["fract_ion-3"].rolling(window=window).std()

    t1_1 = im1["Ti_ion-1"].rolling(window=window).mean()
    t2_1 = im1["Ti_ion-2"].rolling(window=window).mean()
    #t3_1 = im1["Ti_ion-3"].rolling(window=window).mean()
    std_t1_1 = im1["Ti_ion-1"].rolling(window=window).std()
    std_t2_1 = im1["Ti_ion-2"].rolling(window=window).std()
    #std_t3_1 = im1["Ti_ion-3"].rolling(window=window).std()

    v1_1 = im1["Va_ion-1"].rolling(window=window).mean()
    v2_1 = im1["Va_ion-2"].rolling(window=window).mean()
    #v3_1 = im1["Va_ion-3"].rolling(window=window).mean()
    std_v1_1 = im1["Va_ion-1"].rolling(window=window).std()
    std_v2_1 = im1["Va_ion-2"].rolling(window=window).std()
    #std_v3_1 = im1["Va_ion-3"].rolling(window=window).std()


    ax1[1,0].plot(r, f1_1,color='orangered',linewidth=2, linestyle='-', label='H')
    ax1[1,0].fill_between(r,f1_1-std_f1_1,f1_1+std_f1_1,
                color='coral',        # fill color
                alpha=0.2,           # transparency
                edgecolor='tomato',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    #ax1[0].plot(r, lp["fract_ion-2"],color='mediumblue',linewidth=2, linestyle='-', label='Cold H')
    ax1[1,0].plot(r, f2_1,color='MediumSeaGreen',linewidth=2, linestyle='-', label='He')
    ax1[1,0].fill_between(r,f2_1-std_f2_1,f2_1+std_f2_1,
                color='MediumAquamarine',        # fill color
                alpha=0.2,           # transparency
                edgecolor='MediumSeaGreen',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )

    ax1[1,0].set_title('Fractions')
    ax1[1,0].set_xlabel(r"$x (mm)$")
    ax1[1,0].set_xlim([3.55, 4.75])
    ax1[1,0].invert_xaxis()
    ax1[1,0].set_ylabel(r"$Fraction$")
    ax1[1,0].grid(True)
    ax1[1,0].legend()

    ax1[0,1].plot(r, t1_1,color='orangered',linewidth=2,linestyle='-',label='H')
    ax1[0,1].fill_between(r,t1_1-std_t1_1,t1_1+std_t1_1,
                color='coral',        # fill color
                alpha=0.2,           # transparency
                edgecolor='tomato',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    #ax1[1].plot(r, lp["Ti_ion-2"],color='mediumblue',linewidth=2,linestyle='-')
    ax1[0,1].plot(r, t2_1,color='MediumSeaGreen',linewidth=2,linestyle='-',label='He')
    ax1[0,1].fill_between(r,t2_1-std_t2_1,t2_1+std_t2_1,
                color='MediumAquamarine',        # fill color
                alpha=0.2,           # transparency
                edgecolor='MediumSeaGreen',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    ax1[0,1].set_title('Temperatures')
    ax1[0,1].set_xlabel(r"$x (mm)$")
    ax1[0,1].set_xlim([3.55, 4.75])
    ax1[0,1].set_ylim([0.2, 0.85])
    ax1[0,1].invert_xaxis()
    ax1[0,1].set_ylabel(r"$T_{i} (keV)$")
    ax1[0,1].grid(True)
    ax1[0,1].legend()

    ax1[1,1].plot(r, v1_1,color='orangered',linewidth=2,linestyle='-',label='H')
    ax1[1,1].fill_between(r,v1_1-std_v1_1,v1_1+std_v1_1,
                color='coral',        # fill color
                alpha=0.2,           # transparency
                edgecolor='tomato',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    #ax1[2].plot(r, lp["Va_ion-2"],color='mediumblue',linewidth=2,linestyle='-')
    ax1[1,1].plot(r, v2_1,color='MediumSeaGreen',linewidth=2,linestyle='-',label='He')
    ax1[1,1].fill_between(r,v2_1-std_v2_1,v2_1+std_v2_1,
                color='MediumAquamarine',        # fill color
                alpha=0.2,           # transparency
                edgecolor='MediumSeaGreen',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    #ax1[2].plot(r, bv,color='k',linewidth=2,linestyle='--')
    ax1[1,1].set_title('Velocities')
    ax1[1,1].set_xlabel(r"$x (mm)$")
    ax1[1,1].set_xlim([3.55, 4.75])
    ax1[1,1].invert_xaxis()
    ax1[1,1].set_ylabel(r"$V_{i} (10^{6} cm/s)$")
    ax1[1,1].grid(True)
    ax1[1,1].legend()