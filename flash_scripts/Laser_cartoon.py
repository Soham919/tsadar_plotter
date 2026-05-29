import numpy as np
import matplotlib.pyplot as plt

t = np.array([0,0.177,0.634,0.818])
p = np.array([0, 3, 3, 0])

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

fig, ax = plt.subplots(figsize=(8,4))
ax.plot(t, p, color='red', linewidth=2, label=r'Drive Beam : 3$\omega$')
ax.set_xlabel(r"$t (ns)$")
ax.set_ylabel(r"$P (TW)$")
ax.set_title(r'3TW, 600 ps flat Laser Pulse profile', pad=15)
ax.set_xlim([0, 1])
ax.set_ylim([0, 4])
ax.grid(True)
plt.legend()
plt.show()

