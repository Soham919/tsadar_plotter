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
shots = ["92534"]  # change this for the shots you want to plot
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
    #fp = baseDir/".."/"Kinshock"/"Kinshock-19A"/"Data"/shot
    fp = Path("//profiles/Users$/sban/Documents/Kinshock/Kinshock-19A/Data")
    fp = fp/shot
    # --------- Mac ------------
    # fp = baseDir/"Kinshock-19A"/"data"/"117828"

    # saveDir = baseDir/"Kinshock"/"Kinshock-26A"/"Data"/shot/"plots"
    # if not saveDir.exists():
    #     saveDir.mkdir(parents = True, exist_ok = True)


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

# ============================================================
# FIGURE SETTINGS
# ============================================================

FIGSIZE = (5, 5)

FONT_SIZE = 16
LABEL_SIZE = 18
TICK_SIZE = 16

AXIS_WIDTH = 1.2
TICK_WIDTH = 1.2
TICK_LENGTH = 5

PROBE_COLOR = "white"
PROBE_LINEWIDTH = 2.0

plt.rcParams.update({
    "font.family": "serif",
    "font.size": FONT_SIZE,
    "mathtext.fontset": "cm",

    "axes.labelsize": LABEL_SIZE,
    "axes.linewidth": AXIS_WIDTH,

    "xtick.labelsize": TICK_SIZE,
    "ytick.labelsize": TICK_SIZE,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.width": TICK_WIDTH,
    "ytick.major.width": TICK_WIDTH,
    "xtick.major.size": TICK_LENGTH,
    "ytick.major.size": TICK_LENGTH,
})

# ============================================================
# Plot
# ============================================================

fig, ax = plt.subplots(figsize=FIGSIZE)

# Convert spatial calibration to mm
x_plot = shot_data.loc[int(shots[0]), "pointing"] - x_epw[0, :] / 1e3

pcm = ax.pcolormesh(
    x_plot,
    y_epw[0, :],
    epw[0, :, :],
    cmap="turbo",
    shading="auto",
    rasterized=True
)

# Probe wavelength
ax.axhline(
    526.5,
    color=PROBE_COLOR,
    linestyle="--",
    linewidth=PROBE_LINEWIDTH
)

# ============================================================
# Axes
# ============================================================

ax.set_xlim(2.25, 4.0)

# ax.set_xticks([
#     2.25, 2.50, 2.75, 3.00,
#     3.25, 3.50, 3.75, 4.00
# ])

ax.set_xlabel(r"$x$ (mm)")
ax.set_xticks([2.50, 3.00, 3.50, 4.00])   
ax.set_ylabel(r"$\lambda$ (nm)")

# No title -- use the LaTeX caption instead

# ============================================================
# Margins
# ============================================================

fig.subplots_adjust(
    left=0.18,
    right=0.96,
    bottom=0.16,
    top=0.97
)

# ============================================================
# Save
# ============================================================

fig.savefig("EPW_spectrum.pdf")

plt.show()