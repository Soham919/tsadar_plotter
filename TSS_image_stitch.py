import h5py
import numpy as np
import matplotlib.pyplot as plt
from h5_helpers.h5_helper import *
from pathlib import Path
from tsadar_plots.TS_aux import shot_data
from tsadar_plots.TS_aux import get_calibration

baseDir = Path().parent.resolve()  # get the parent directory of the current file

#shot = "117827"  # change this for the shot you want to plot
ts_type = "imaging" # change this for the type of TS data you want to plot (imaging or time resolved)
shots = ["117830"]  # change this for the shots you want to plot
epw = np.zeros((len(shots), 1024, 1024))
iaw = np.zeros((len(shots), 1024, 1024))

# ---- variables for axes --- #
x_epw = np.zeros((len(shots), 1024))
x_iaw = np.zeros((len(shots), 1024))
y_epw = np.zeros((len(shots), 1024))
y_iaw = np.zeros((len(shots), 1024))

i = 0
for shot in shots:
    # --------- Windows ------------ #
    fp = baseDir/".."/"Kinshock"/"Kinshock-26A"/"Data"/shot
    # --------- Mac ------------
    # fp = baseDir/"Kinshock-19A"/"data"/"117828"

    saveDir = baseDir/"Kinshock"/"Kinshock-26A"/"Data"/shot/"plots"
    if not saveDir.exists():
        saveDir.mkdir(parents = True, exist_ok = True)


    f_iaw = f"IAW_CCD-s{shot}.h5"  # switch names for ROSS and CCD since they are named differently
    f_epw = f"EPW_CCD-s{shot}.h5"  # switch names for ROSS and CCD since they are named differently
    # f_tbd = f"P9TBD_CCD-s{shot}_ccd.h5"
    # f_tpdi = f"TPDI-s{shot}.h5"
    # f_xrphc = f"XRPHC-CID-xphc_{shot}_h8_cid.h5"
    # fp = fp 

    #h5_helper.h5_show(fp/ f_xrphc)

    with h5py.File(fp/ f_iaw,'r') as f:
        d = f["Streak_array"][:]
        d = d[0,:,:] #- d[1,:,:]  # subtract the two streaks to get the IAW signal
        # ------------- Readjust for ccd imgs -------------- #
        d = np.rot90(d, 1)  # rotate the image to be vertical
        d = np.fliplr(d)  # flip the image
        #d = np.flipud(d)  # flip the image 
        iaw[i,:,:] = d

    with h5py.File(fp/ f_epw,'r') as f:
        d = f["Streak_array"][:]
        d = d[0,:,:] #- d[1,:,:]  # subtract the two streaks to get the EPW signal
        # ------------- Readjust for ccd imgs -------------- #
        d = np.rot90(d, 1)  # rotate the image to be vertical
        d = np.fliplr(d)  # flip the image
        d = np.flipud(d)  # flip the image 
        epw[i,:,:] = d
    
   
    
    # --- get calibrated axes for plotting --- #
    x_epw[i,:], x_iaw[i,:], y_epw[i,:], y_iaw[i,:] = get_calibration(shot, ts_type)

    i = i + 1

fig, ax = plt.subplots(1,len(shots),figsize=(1.25*5*len(shots),5), constrained_layout=True)

plt.rcParams.update({
    #    "font.family": "serif",             # pick font family
    #    "font.serif": ["Arial"],  # pick specific font
        "font.size": 60,                     # set size
        "axes.labelsize": 60,    # font size for axis labels
        "axes.titlesize": 60,     # font size for titles
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

for i in range(len(shots)):
    x_iaw[i,:] = shot_data.loc[int(shots[i]),"pointing"] - x_iaw[i,:]/10**3

    if len(shots) == 1:
        ax = [ax]  # make ax a list if there's only one subplot
    ax[i].pcolormesh(x_iaw[i,:], y_iaw[i,:], iaw[i,:,:], cmap='turbo', shading='auto')   # change to iaw/epw if you want to plot the epws
    ax[i].hlines(
    526.5,
    xmin=x_iaw[i,:].min(),
    xmax=x_iaw[i,:].max(),
    color='white',
    linestyle='--',
    linewidth=2
    )
    ax[i].set_xlim([7.2, 8.6])
    ax[i].set_xlabel("x (mm)")
    ax[i].set_ylabel(r"$\lambda (nm)$")

    #ax[i].set_title(f"Shot {shots[i]}")

# Remove some tics for a stitched picture, not needed if plotting just one shot
#ax[1].set_yticks([])      # removes tick marks
#ax[1].set_ylabel("")      # removes label

#ax[2].set_yticks([])      # removes tick marks
#ax[2].set_ylabel("")      # removes label

#plt.tight_layout()
plt.show()