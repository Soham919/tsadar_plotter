import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import h5_helper as h5h

baseDir = Path().resolve().parent
fp = baseDir/"Kinshock"/"Kinshock-26A"/"data"/"Alignment"
fn = "IAWCCD_2w_532cw_100ph.hdf"

h5h.h5_show(fp/fn)
img = h5h.h5_import(fp/fn, "Streak_array")

fig, ax  = plt.subplots()
im = ax.imshow(img[0,:,:], cmap='jet')  # adjust vmin and vmax for better contrast
plt.colorbar(im, label='counts')
plt.show()

