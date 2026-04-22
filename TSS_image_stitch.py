import h5py
import numpy as np
import matplotlib.pyplot as plt
import h5_helper
from pathlib import Path
from TS_aux import shot_data
from TS_aux import get_calibration

# baseDir = Path(__file__).resolve().parent

# #shot = "117827"  # change this for the shot you want to plot
# ts_type = "imaging" # change this for the type of TS data you want to plot (imaging or time resolved)
# shots = ["92535", "92536", "92537"]
# epw = np.zeros((len(shots), 1024, 1024))
# iaw = np.zeros((len(shots), 1024, 1024))

# # ---- variables for axes --- #
# x_epw = np.zeros((len(shots), 1024))
# x_iaw = np.zeros((len(shots), 1024))
# y_epw = np.zeros((len(shots), 1024))
# y_iaw = np.zeros((len(shots), 1024))

# i = 0
# for shot in shots:
#     # --------- Windows ------------ #
#     fp = baseDir/"Kinshock"/"Kinshock-26A"/"data"/shot
#     # --------- Mac ------------
#     # fp = baseDir/"Kinshock-26A"/"data"/"117828"

#     saveDir = baseDir/"Kinshock"/"Kinshock-26A"/"data"/shot/"plots"
#     if not saveDir.exists():
#         saveDir.mkdir(parents = True, exist_ok = True)


#     f_iaw = f"IAW_CCD-s{shot}.h5"  # switch names for ROSS and CCD since they are named differently
#     f_epw = f"EPW_CCD-s{shot}.h5"  # switch names for ROSS and CCD since they are named differently
#     # f_tbd = f"P9TBD_CCD-s{shot}_ccd.h5"
#     # f_tpdi = f"TPDI-s{shot}.h5"
#     # f_xrphc = f"XRPHC-CID-xphc_{shot}_h8_cid.h5"
#     # fp = fp 

#     #h5_helper.h5_show(fp/ f_xrphc)

#     with h5py.File(fp/ f_iaw,'r') as f:
#         d = f["Streak_array"][:]
#         d = d[0,:,:] - d[1,:,:]  # subtract the two streaks to get the IAW signal
#         iaw[i,:,:] = d

#     with h5py.File(fp/ f_epw,'r') as f:
#         d = f["Streak_array"][:]
#         d = d[0,:,:] - d[1,:,:]  # subtract the two streaks to get the EPW signal
#         epw[i,:,:] = d
    
#     # --- get calibrated axes for plotting --- #
#     x_epw[i,:], x_iaw[i,:], y_epw[i,:], y_iaw[i,:] = get_calibration(shot, ts_type)

#     i = i + 1

# # ------------- Readjust for ccd imgs -------------- #
# # epw = np.rot90(epw[0,:,:], 1)  # rotate the image to be vertical
# # iaw = np.rot90(iaw[0,:,:], 1)  # rotate the image to be vertical
# # iaw = np.fliplr(iaw)  # flip the image
# # epw = np.fliplr(epw)  # flip the image
# # iaw = np.flipud(iaw)  # flip the image 

# fig, ax = plt.subplots(1,3,figsize=(30,10))

# for i in range(len(shots)):
#     ax[i].imshow(x_iaw[i,:], y_iaw[i,:], iaw[i,:,:], aspect='auto')
#     ax[i].hlines(526.5, color='white', linestyle='--', linewidth=1)
#     ax[i].set_title(f"IAW - Shot {shots[i]}")
#     ax[i].set_xlabel("x (mm)")
#     ax[i].set_ylabel(r"$\lambda (nm)$")