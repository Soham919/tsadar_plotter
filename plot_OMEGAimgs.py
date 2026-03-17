import h5py
import numpy
import matplotlib.pyplot as plt
import h5_helper
from pathlib import Path

baseDir = Path().resolve().parent
fp = baseDir/"Kinshock-26A"/"data"/"117828"
file = "IAW-s117828.h5"
fp = fp / file

h5_helper.h5_show(fp)

with h5py.File(fp,'r') as f:
    img = f["Streak_array"][:]

fig, ax  = plt.subplots()
ax = plt.imshow(img[0,:,:])
plt.show()
