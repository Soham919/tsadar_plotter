import pandas as pd
import xarray as xr
from dataclasses import dataclass
from datetime import date
import numpy as np
import matplotlib.pyplot as plt
import cv2
import os 
from pathlib import Path
import re


shot_data = pd.DataFrame({
    "shot" : [92522,92524,92525,92527,92528,92530,92531,92532,92533,92534,92535,92536,92537,92538,117830,117831,117832,117837,117838,117839,117827,117828],
    "pointing" : [4,2,4,6,8,5.2,4,4,5,3,6,5,4,4,8,8,5.5,5.5,5.5,5.5,8,8],
    "timing" : [4,2,4,6,8,6,4,4,4,4,6,6,6,4,9,9,7.6,7.3,6,5,0,0]
})

shot_data.set_index("shot", inplace=True)

def load_file(data, file_tag, file_path):
    """
    Load multiple CSV files for a given shot number into a nested dictionary.
    """
    try:
        s = file_tag.split('.')
        if s[-1]=="csv":
            data[file_tag] = pd.read_csv(file_path)
            print(f"Loaded {file_tag}")

        if s[-1]=="nc":
            data[file_tag] = xr.open_dataset(file_path)
            print(data[file_tag])
            
    except FileNotFoundError:
        print(f"⚠ File not found: {file_path}")

def select_shots(folder_path):
    """
    Select specific shots from the data dictionary.
    """
    ls = [f for f in os.listdir(folder_path) if os.path.isdir(Path(folder_path)/f)]
    data = {}

    print("Stuff in folder :\n")
    for i in range(len(ls)):
        print(f"[{i:g}] {ls[i]}\n")
        
    while True:
        n = int(input("Which shot would you like to see fits for?\n"))
        fp2 = os.listdir(Path(folder_path)/ls[n])
        if str(ls[n]) not in data:
            data[str(ls[n])] = {}
        lscsv = [f for f in fp2]
        print(f"Files for shot {ls[n]} :\n")
        for i in range(len(lscsv)):
            print(f"[{i}] {lscsv[i]}")
        
        fn = input("Which files do you want to import?\n")
        files  = re.split(r"[ ,]+", fn)
        files = [int(x) for x in files]
        for i in files:
            fp = Path(folder_path)/ls[n]/lscsv[i]
            if fp.is_dir():
                if str(lscsv[i]) not in data[str(ls[n])]:
                    data[str(ls[n])][str(lscsv[i])] = {}
                files2 = os.listdir(fp)
                for f in files2:
                    load_file(data[str(ls[n])][str(lscsv[i])],f,Path(fp)/f)
            else:
                load_file(data[str(ls[n])],lscsv[i],fp)
                
        flag = input("Load more shots?\n")
        if flag in ("n","no","nah","","nope"):
            break

    return data

def plot_lineouts(data):
    #### PLOT LINEOUTS AND FITS ####
    bin_file = data['117830-6']['ion_fit_and_data.nc']
    bin_file2 = data['117830-10']['ion_fit_and_data.nc']
    bin_file3 = data['117830-11']['ion_fit_and_data.nc']
    #e_file = data[ls[n]]['92536_13_ele_fit_and_data.nc']
    #e_file2 = data[ls[n]]['92536_15_ele_fit_and_data.nc']
    img = bin_file["data"].values
    fit = bin_file["fit"].values
    #e_img = e_file["data"].values
    #e_fit = e_file["data"].values

    img2 = bin_file2["data"].values
    fit2 = bin_file2["fit"].values
    #e2_img = e_file2["data"].values
    #e2_fit = e_file2["data"].values

    img3 = bin_file3["data"].values
    fit3 = bin_file3["fit"].values

    r = bin_file[r"Radius (\mum)"].values
    l = bin_file["Wavelength"].values

    r2 = bin_file2[r"Radius (\mum)"].values
    l2 = bin_file2["Wavelength"].values

    r3 = bin_file3[r"Radius (\mum)"].values
    l3 = bin_file3["Wavelength"].values

    x = int(input("Which radius do you want to plot lineout for?\n"))
    idx = (np.abs(r-x)).argmin()
    idx2 = (np.abs(r2-x)).argmin()
    idx3 = (np.abs(r3-x)).argmin()
    fig, ax = plt.subplots(figsize=(8,5))
    # ax.plot(l2,img2[idx2,:],lw = 2,color='orange',label=f'data({int(r[idx])})')
    # ax.plot(l,fit[idx,:],lw = 2,color='orangered',linestyle='--', label=f'2 species fit({int(r[idx])})')
    # ax.plot(l2,fit2[idx2,:],lw=2,color='mediumblue',linestyle='--', label=f'3 species fit({int(r2[idx2])})')

    # For the TS img plot
    ax.plot(l2,img2[idx2,:],lw=1, color='orange',label=f' 117830 Thomson Data')
    ax.plot(l,fit[idx,:],lw = 2,color='blue',linestyle='--', label=f'Z=2')
    ax.plot(l2,fit2[idx2,:],lw=2,color='red',linestyle='--', label=f'Z=1')
    ax.plot(l3,fit3[idx3,:],lw=2,color='green',linestyle='--', label=f'Z=1+2')

    # ax.plot(l,fit[idx+1,:],color='blue',linestyle='--', label=f'fit({r[idx+1]})')
    # ax.plot(l,fit[idx-1,:],color='green',linestyle='--', label=f'fit({r[idx-1]})')
    ax.set_xlabel(r"$\lambda (nm)$")
    ax.set_ylabel(r"$Amplitude (a.u.)$")
    ax.set_xlim([524.0, 528.0])
    ax.set_title(f"Lineout at r={r[idx]:.1f} um")
    #ax.set_ylim([3.55, 4.75])

    # plt.rcParams.update({
    #    'font.family': 'serif',
    #    'font.serif': ['Times New Roman'],
    #     'font.size': 10,
    #     'axes.labelsize': 11,
    #     'axes.titlesize': 11,
    #     'xtick.labelsize': 9,
    #     'ytick.labelsize': 9,
    #     'legend.fontsize': 9,
    #  })

    ax.tick_params(direction='in', length=3.5, width=0.8, colors='black')
    ax.xaxis.set_ticks_position('both')
    ax.yaxis.set_ticks_position('both')

    plt.legend()


def get_calibration(shotNum,tstype):
    '''
    Returns the calibration parameters for the EPW and IAW spectrometers based on the shot number and type of spectrometer.
    Parameters:
    shotNum (int): The shot number for which to retrieve the calibration parameters.
    tstype (str): The type of TS imaging for the shot, either "imaging" or "time resolved", current;y only works for "imaging".
    Returns:
    tuple: A tuple containing the x and y axis of the EPW and IAW spectrometers, respectively.
    '''
    shotNum = int(shotNum)
    CCDsize = (1024, 1024)
    if shotNum in [92522,92525,92531,92532,92534,92537,92538]:
        EPWDisp = 0.27093
        IAWDisp = 0.0057
        EPWoff = 384.80361  # only this shot seems to have shifted
        IAWoff = 523.74

        #stddev["spect_stddev_ion"] = 0.028  # needs to be checked
        #stddev["spect_stddev_ele"] = 1.4365  # needs to be checked

        magI = 2.87 #3.8  #2.87  # um / px
        magE = 5.13  # um / px

        EPWtcc = 565 #1024 - 456.1  # 562;
        IAWtcc = 625 #519#1024 - 519  # 469;
            
    elif shotNum in [92527,92535]:
        EPWDisp = 0.27093
        IAWDisp = 0.0057
        EPWoff = 384.53268  
        IAWoff = 523.6842

        #stddev["spect_stddev_ion"] = 0.028  # needs to be checked
        #stddev["spect_stddev_ele"] = 1.4365  # needs to be checked

        magI = 2.87  # um / px
        magE = 5.13  # um / px

        EPWtcc = 565 #1024 - 456.1  # 562;
        IAWtcc = 632 #519#1024 - 519  # 469;

    elif shotNum in [92528,92533,92536]:
        EPWDisp = 0.27093
        IAWDisp = 0.0057
        EPWoff = 385.80361 #379.38500  
        IAWoff = 523.6899

        #stddev["spect_stddev_ion"] = 0.028  # needs to be checked
        #stddev["spect_stddev_ele"] = 1.4365  # needs to be checked

        magI = 2.87  # um / px
        magE = 5.13  # um / px

        EPWtcc = 563 #1024 - 456.1  # 562;
        IAWtcc = 630 #519#1024 - 519  # 469;

    elif 117830 <= shotNum <= 117839:
        EPWDisp = 0.276
        IAWDisp = 0.00437
        EPWoff = 375.528  # needs to be checked
        IAWoff = 524.189

        #stddev["spect_stddev_ion"] = 0.028  # needs to be checked
        #stddev["spect_stddev_ele"] = 1.4365  # needs to be checked

        magI = 2.89  # um / px times strech factor accounting for tilt in view
        magE = 5.13 # um / px times strech factor accounting for tilt in view

        EPWtcc = 528  # 562;
        IAWtcc = 475  # 469;

    ## Apply calibrations
    axisx = np.arange(1, CCDsize[1] + 1)
    axisxE = (axisx) * magE  # ps,um
    axisxI = (axisx) * magI  # ps,um
    axisy = np.arange(1, CCDsize[0] + 1)
    axisyE = axisy * EPWDisp + EPWoff  # (nm)
    axisyI = axisy * IAWDisp + IAWoff  # (nm)

    if tstype == "imaging":
        axisxE = axisxE - EPWtcc * magE
        axisxI = axisxI - IAWtcc * magI

    return axisxE, axisxI, axisyE, axisyI