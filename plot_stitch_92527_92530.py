import numpy as np
import matplotlib.pyplot as plt
import cv2
import os 
from pathlib import Path
import re
import pandas as pd
from TS_aux import load_file
from TS_aux import shot_data

def plot_stitch_92527_92530(data):
    ### PLOT STITCHED H2 SHOCKS ###
    im1 = data['92527']['92527_19_nersc.csv']
    im2 = data['92530']['92530_7_nersc.csv']
    #im3 = data['92537']['92537_5_nersc.csv']
    magI = 2.87
    r1 = shot_data.loc[92527,"pointing"] - im1["Radius (\\mum)"]/10**3 + ((620-519)*magI/10**3) # convert to mm and shift 
    r2 = shot_data.loc[92530,"pointing"] - im2["Radius (\\mum)"]/10**3 + ((632-519)*magI/10**3)
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

    fig1, ax1  = plt.subplots(2,2,figsize=(16,13))
    # Error parameters
    window = 3  # rolling window
    f1_1 = im1["fract_ion-1"].rolling(window=window).mean()
    f2_1 = im1["fract_ion-2"].rolling(window=window).mean()
    f3_1 = im1["fract_ion-3"].rolling(window=window).mean()
    std_f1_1 = im1["fract_ion-1"].rolling(window=window).std()
    std_f2_1 = im1["fract_ion-2"].rolling(window=window).std()
    std_f3_1 = im1["fract_ion-3"].rolling(window=window).std()

    f1_2 = im2["fract_ion-1"].rolling(window=window).mean()
    f2_2 = im2["fract_ion-2"].rolling(window=window).mean()
    f3_2 = im2["fract_ion-3"].rolling(window=window).mean()
    std_f1_2 = im2["fract_ion-1"].rolling(window=window).std()
    std_f2_2 = im2["fract_ion-2"].rolling(window=window).std()
    std_f3_2 = im2["fract_ion-3"].rolling(window=window).std()

    t1_1 = im1["Ti_ion-1"].rolling(window=window).mean()
    t2_1 = im1["Ti_ion-2"].rolling(window=window).mean()
    t3_1 = im1["Ti_ion-3"].rolling(window=window).mean()
    std_t1_1 = im1["Ti_ion-1"].rolling(window=window).std()
    std_t2_1 = im1["Ti_ion-2"].rolling(window=window).std()
    std_t3_1 = im1["Ti_ion-3"].rolling(window=window).std()

    t1_2 = im2["Ti_ion-1"].rolling(window=window).mean()
    t2_2 = im2["Ti_ion-2"].rolling(window=window).mean()
    t3_2 = im2["Ti_ion-3"].rolling(window=window).mean()
    std_t1_2 = im2["Ti_ion-1"].rolling(window=window).std()
    std_t2_2 = im2["Ti_ion-2"].rolling(window=window).std()
    std_t3_2 = im2["Ti_ion-3"].rolling(window=window).std()

    v1_1 = im1["Va_ion-1"].rolling(window=window).mean()
    v2_1 = im1["Va_ion-2"].rolling(window=window).mean()
    v3_1 = im1["Va_ion-3"].rolling(window=window).mean()
    std_v1_1 = im1["Va_ion-1"].rolling(window=window).std()
    std_v2_1 = im1["Va_ion-2"].rolling(window=window).std()
    std_v3_1 = im1["Va_ion-3"].rolling(window=window).std()

    v1_2 = im2["Va_ion-1"].rolling(window=window).mean()
    v2_2 = im2["Va_ion-2"].rolling(window=window).mean()
    v3_2 = im2["Va_ion-3"].rolling(window=window).mean()
    std_v1_2 = im2["Va_ion-1"].rolling(window=window).std()
    std_v2_2 = im2["Va_ion-2"].rolling(window=window).std()
    std_v3_2 = im2["Va_ion-3"].rolling(window=window).std()



    ax1[1,0].plot(r1, f1_1,color='orangered',linewidth=2, linestyle='-', label='hot H')
    ax1[1,0].plot(r1, f2_1,color='mediumblue',linewidth=2, linestyle='-', label='cold H')
    ax1[1,0].plot(r1, f3_1,color='MediumSeaGreen',linewidth=2, linestyle='-', label='N')
    ax1[1,0].fill_between(r1,f1_1-std_f1_1,f1_1+std_f1_1,
                color='coral',        # fill color
                alpha=0.2,           # transparency
                edgecolor='tomato',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    ax1[1,0].fill_between(r1,f2_1-std_f2_1,f2_1+std_f2_1,
                color='CornflowerBlue',        # fill color
                alpha=0.2,           # transparency
                edgecolor='mediumblue',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    ax1[1,0].fill_between(r1,f3_1-std_f3_1,f3_1+std_f3_1,
                color='MediumAquamarine',        # fill color
                alpha=0.2,           # transparency
                edgecolor='MediumSeaGreen',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )


    ax1[1,0].plot(r2, f1_2,color='orangered',linewidth=2, linestyle='--')
    ax1[1,0].plot(r2, f2_2,color='mediumblue',linewidth=2, linestyle='--')
    ax1[1,0].plot(r2, f3_2,color='MediumSeaGreen',linewidth=2, linestyle='--')
    ax1[1,0].fill_between(r2,f1_2-std_f1_2,f1_2+std_f1_2,
                color='coral',        # fill color
                alpha=0.2,           # transparency
                edgecolor='orangered',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    ax1[1,0].fill_between(r2,f2_2-std_f2_2,f2_2+std_f2_2,
                color='CornflowerBlue',        # fill color
                alpha=0.2,           # transparency
                edgecolor='mediumblue',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    ax1[1,0].fill_between(r2,f3_2-std_f3_2,f3_2+std_f3_2,
                color='MediumAquamarine',        # fill color
                alpha=0.2,           # transparency
                edgecolor='MediumSeaGreen',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    #ax1[0].plot(r3, im3["fract_ion-1"],color='orangered',linewidth=2, linestyle=':')
    #ax1[0].plot(r3, im3["fract_ion-2"],color='mediumblue',linewidth=2, linestyle=':')
    #ax1[0].plot(r3, im3['fract_ion-2'],color='MediumSeaGreen',linewidth=2, linestyle=':')

    ax1[1,0].set_title('Fractions')
    ax1[1,0].set_xlabel(r"$x (mm)$")
    ax1[1,0].set_xlim([4.6,6.6])
    ax1[1,0].set_xticks([4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6])
    #ax1[1,0].tick_params(labelbottom=False)
    ax1[1,0].set_ylabel(r"$Fraction$")
    ax1[1,0].legend()
    ax1[1,0].grid(True)


    ax1[0,1].plot(r1, t1_1,color='orangered',linewidth=2, linestyle='-')
    ax1[0,1].plot(r1, t2_1,color='mediumblue',linewidth=2, linestyle='-')
    ax1[0,1].plot(r1, t3_1,color='mediumseagreen',linewidth=2, linestyle='-')
    ax1[0,1].fill_between(r1,t1_1-std_t1_1,t1_1+std_t1_1,
                color='coral',        # fill color
                alpha=0.2,           # transparency
                edgecolor='tomato',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    ax1[0,1].fill_between(r1,t2_1-std_t2_1,t2_1+std_t2_1,
                color='CornflowerBlue',        # fill color
                alpha=0.2,           # transparency
                edgecolor='mediumblue',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    ax1[0,1].fill_between(r1,t3_1-std_t3_1,t3_1+std_t3_1,
                color='MediumAquamarine',        # fill color
                alpha=0.2,           # transparency
                edgecolor='MediumSeaGreen',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )


    ax1[0,1].plot(r2, t1_2,color='orangered',linewidth=2, linestyle='--')
    ax1[0,1].plot(r2, t2_2,color='mediumblue',linewidth=2, linestyle='--')
    ax1[0,1].plot(r2, t3_2,color='MediumSeaGreen',linewidth=2, linestyle='--')
    ax1[0,1].fill_between(r2,t1_2-std_t1_2,t1_2+std_t1_2,
                color='coral',        # fill color
                alpha=0.2,           # transparency
                edgecolor='tomato',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    ax1[0,1].fill_between(r2,t2_2-std_t2_2,t2_2+std_t2_2,
                color='CornflowerBlue',        # fill color
                alpha=0.2,           # transparency
                edgecolor='mediumblue',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    ax1[0,1].fill_between(r2,t3_2-std_t3_2,t3_2+std_t3_2,
                color='MediumAquamarine',        # fill color
                alpha=0.2,           # transparency
                edgecolor='MediumSeaGreen',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )

    #ax1[1].plot(r3, im3["Ti_ion-1"],color='orangered',linewidth=2, linestyle=':')
    #ax1[1].plot(r3, im3["Ti_ion-2"],color='mediumblue',linewidth=2, linestyle=':')
    #ax1[1].plot(r3, im3["Ti_ion-2"],color='MediumSeaGreen',linewidth=2, linestyle=':')
    ax1[0,1].legend()
    ax1[0,1].grid(True)
    ax1[0,1].set_title('Temperatures')
    ax1[0,1].set_xlabel(r"$x (mm)$")
    ax1[0,1].set_xlim([4.6,6.6])
    ax1[0,1].set_xticks([4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6])
    #ax1[0,1].tick_params(labelbottom=False)
    ax1[0,1].set_ylabel(r"$T_{i} (KeV)$")


    ax1[1,1].plot(r1, v1_1,color='orangered',linewidth=2, linestyle='-')
    ax1[1,1].plot(r1, v2_1,color='mediumblue',linewidth=2, linestyle='-',label='TS IAW 2')
    ax1[1,1].plot(r1, v3_1,color='MediumSeaGreen',linewidth=2, linestyle='-')
    ax1[1,1].fill_between(r1,v1_1-std_v1_1,v1_1+std_v1_1,
                color='coral',        # fill color
                alpha=0.2,           # transparency
                edgecolor='tomato',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    ax1[1,1].fill_between(r1,v2_1-std_v2_1,v2_1+std_v2_1,
                color='CornflowerBlue',        # fill color
                alpha=0.2,           # transparency
                edgecolor='mediumblue',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    ax1[1,1].fill_between(r1,v3_1-std_v3_1,v3_1+std_v3_1,
                color='MediumAquamarine',        # fill color
                alpha=0.2,           # transparency
                edgecolor='MediumSeaGreen',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )



    ax1[1,1].plot(r2, v1_2,color='orangered',linewidth=2, linestyle='--')
    ax1[1,1].plot(r2, v2_2,color='mediumblue',linewidth=2, linestyle='--',label='TS IAW 1')
    ax1[1,1].plot(r2, v3_2,color='MediumSeaGreen',linewidth=2, linestyle='--')
    ax1[1,1].fill_between(r2,v1_2-std_v1_2,v1_2+std_v1_2,
                color='coral',        # fill color
                alpha=0.2,           # transparency
                edgecolor='tomato',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    ax1[1,1].fill_between(r2,v2_2-std_v2_2,v2_2+std_v2_2,
                color='CornflowerBlue',        # fill color
                alpha=0.2,           # transparency
                edgecolor='mediumblue',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )
    ax1[1,1].fill_between(r2,v3_2-std_v3_2,v3_2+std_v3_2,
                color='MediumAquamarine',        # fill color
                alpha=0.2,           # transparency
                edgecolor='MediumSeaGreen',    # outline color
                linewidth=1.5,       # outline width
                linestyle='-',      # outline style (optional)
    )


    #ax1[2].plot(r3, im3["Va_ion-1"],color='orangered',linewidth=2, linestyle=':')
    #ax1[2].plot(r3, im3["Va_ion-2"],color='mediumblue',linewidth=2, linestyle=':')
    #ax1[2].plot(r3, im3["Va_ion-2"],color='MediumSeaGreen',linewidth=2, linestyle=':')
    ax1[1,1].set_title('Velocities')
    ax1[1,1].set_xlabel(r"$x (mm)$")
    ax1[1,1].set_xlim([4.6,6.6])
    ax1[1,1].set_xticks([4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6])
    ax1[1,1].set_ylabel(r"$V_{i} (\times 10^{6} cm/s)$")
    ax1[1,1].legend()
    ax1[1,1].grid(True)