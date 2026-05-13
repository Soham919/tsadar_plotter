import h5py
import numpy as np
import matplotlib.pyplot as plt
import h5_helper
from pathlib import Path

baseDir = Path().resolve().parent
# --------- Mac ------------
# fp = baseDir/"Kinshock-26A"/"data"/"117828"
# file = "IAW-s117828.h5"
# fp = fp / file

# --------- Windows ------------
shot = "92532"  # change this for the shot you want to plot
fp = baseDir/"Kinshock"/"Kinshock-19A"/"data"/shot
saveDir = baseDir/"Kinshock"/"Kinshock-19A"/"data"/shot/"plots"
if not saveDir.exists():
    saveDir.mkdir(parents = True, exist_ok = True)


f_iaw = f"IAW-s{shot}.h5"  # switch names for ROSS and CCD since they are named differently
f_epw = f"EPW-s{shot}.h5"  # switch names for ROSS and CCD since they are named differently
f_tbd = f"P9TBD_CCD-s{shot}_ccd.h5"
f_tpdi = f"TPDI-s{shot}.h5"
f_xrphc = f"XRPHC-CID-xphc_{shot}_h8_cid.h5"
fp = fp 

h5_helper.h5_show(fp/ f_xrphc)

with h5py.File(fp/ f_iaw,'r') as f:
    iaw = f["Streak_array"][:]

with h5py.File(fp/ f_epw,'r') as f:
    epw = f["Streak_array"][:]

#with h5py.File(fp/ f_tbd,'r') as f:
#    tbd = f["Streak_array"][:]

with h5py.File(fp/ f_tpdi,'r') as f:
    tpdi = f["Streak_array"][:]

#with h5py.File(fp/ f_xrphc,'r') as f:
#    xr_bg = f["cid_background"][:]
#    xr_fg = f["cid_foreground"][:]

## ------- Plot the XRPHC image --------
# fig, ax  = plt.subplots()
# xr_fg = xr_fg - xr_bg
# xr_fg = np.rot90(xr_fg, 2)  # rotate the image to be vertical
# im = ax.imshow(xr_fg, cmap='grey')  # adjust vmin and vmax for better contrast
# plt.colorbar(im, label='counts')
# ax.set_xticks([])
# ax.set_yticks([])
# ax.set_title(f"Shot {shot} - XRPHC-CID")
# save_path = saveDir / f"shot_{shot}_xrphc.png"
# plt.savefig(save_path, dpi=300, bbox_inches='tight')
# plt.show()

## ------- Plot the TBD image --------
# ax.clear()
# fig, ax  = plt.subplots()
# im = ax.imshow(tbd[0,:,:], cmap='turbo',vmax = 2500)  # adjust vmin and vmax for better contrast
# plt.colorbar(im, label='counts')
# ax.set_xticks([])
# ax.set_yticks([])
# ax.set_title(f"Shot {shot} - P9TBD_CCD")
# save_path = saveDir / f"shot_{shot}_tbd.png"
# plt.savefig(save_path, dpi=300, bbox_inches='tight')
# fig.canvas.draw()
# plt.show()

## ------- Plot the TPDI image --------
#ax.clear()
fig, ax  = plt.subplots()
im = ax.imshow(tpdi[0,:,:], cmap='turbo')  # adjust vmin and vmax for better contrast
plt.colorbar(im, label='counts')
ax.set_xticks([])
ax.set_yticks([])
ax.set_title(f"Shot {shot} - TPDI")
save_path = saveDir / f"shot_{shot}_tpdi.png"
plt.savefig(save_path, dpi=300, bbox_inches='tight')
fig.canvas.draw()
plt.show()


## ------- Plot the IAW and EPW images --------
ax.clear()
fig, ax  = plt.subplots(1,2, figsize=(12,5))

# ------------- Readjust for ccd imgs -------------- #
epw = np.rot90(epw[0,:,:], 1)  # rotate the image to be vertical
iaw = np.rot90(iaw[0,:,:], 1)  # rotate the image to be vertical
iaw = np.fliplr(iaw)  # flip the image
epw = np.fliplr(epw)  # flip the image
iaw = np.flipud(iaw)  # flip the image 
# -------------------------------------------------- #
# ------------- Don't adjust for ROSS imgs --------- #
#epw = epw[0,:,:]
#iaw = iaw[0,:,:]
# -------------------------------------------------- #
im1 = ax[0].imshow(epw, cmap='turbo',vmax=2500)  # adjust vmin and vmax for better contrast
plt.colorbar(im1)
ax[0].set_title(f"Shot {shot} - EPW")

im2 = ax[1].imshow(iaw, cmap='turbo',vmax=2000)  # adjust vmin and vmax for better contrast
plt.colorbar(im2)
ax[1].set_title(f"Shot {shot} - IAW")

save_path = saveDir / f"shot_{shot}_TSS.png"
plt.savefig(save_path, dpi=300, bbox_inches='tight')
fig.canvas.draw()
plt.show()

