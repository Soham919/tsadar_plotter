import numpy as np
import h5py
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import constants
from h5_import import h5_to_dict

#---------File path---------- #
baseDir = Path().resolve().parent
fp = baseDir/"Kinshock-26A"/"GasJet"
file = "D-GJ-C-232_H2He_1500psi.h5"
fp = fp/file

# ---------Load H5 file---------- #
with h5py.File(fp, 'r') as f:
    full_data = h5_to_dict(f)
    #f.visititems(print_h5)

# Density contours at specific z/r locations
target = [6.5e-3, 7.25e-3, 8e-3, 10e-3, 15e-3]  # mm in meters
arr = full_data['r']
idx = np.zeros(len(target), dtype=int)
for i in range(len(target)):
    idx[i] = np.abs(arr - target[i]).argmin()

# ------ Data to visualize -------- #
q = full_data['total-pressure']  # radial velocity in m/s

# ------ Plot -------- #
fig, ax = plt.subplots(1,2, figsize=(12,5))
im = ax[0].pcolormesh(full_data['r']*1000, full_data['z']*1000, q, cmap='viridis', vmin=1e6, vmax = 1e7,shading='auto')
ax[0].plot(np.linspace(6.5,6.5,len(arr)), full_data['z']*1000, color='red', linestyle='--', label='6.5 mm')
ax[0].plot(np.linspace(7.25,7.25,len(arr)), full_data['z']*1000, color='red', linestyle='--', label='7.25 mm')
ax[0].plot(np.linspace(8,8,len(arr)), full_data['z']*1000, color='red', linestyle='--', label='8 mm')
ax[0].plot(np.linspace(10,10,len(arr)), full_data['z']*1000, color='red', linestyle='--', label='10 mm')
ax[0].plot(np.linspace(15,15,len(arr)), full_data['z']*1000, color='red', linestyle='--', label='15 mm')
ax[0].legend()
ax[0].set_xlabel('r (mm)')
ax[0].set_ylabel('z (mm)')
plt.colorbar(im, label='total-pressure (Pa)', ax=ax[0])

for i in idx:
    ax[1].plot(full_data['z']*1000, q[:,i], linestyle='-',label=f'r = {full_data["r"][i]*1000} mm')
ax[1].axvline(x=5, color='grey', linestyle='--', label='z = 5 mm')
ax[1].axvline(x=2.5, color='black', linestyle='--', label='z = 5 mm')
ax[1].axvline(x=7.5, color='black', linestyle='--', label='z = 5 mm')
ax[1].set_xlabel('z (mm)')
ax[1].set_ylabel('Pressure (Pa)')
ax[1].legend()
plt.show()
