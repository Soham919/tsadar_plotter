import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

class LaserEnergy:
    def __init__(self, data):
        self.step = data[:, 0]
        self.time = data[:, 1]
        self.dt = data[:, 2]
        self.ein = data[:, 3]/10**7 # convert from erg to J
        self.eout = data[:, 4]/10**7
        self.dein = data[:, 5]/10**7
        self.deout = data[:, 6]/10**7

        self.eabs = self.ein - self.eout
        self.abs_frac = self.eabs / self.ein

# -------------------------
# Example usage
# -------------------------
local_dir = Path(r"C:\Simulation_data\FLASH\2D_Cylindrical\debug")
run = ["2DCartesian_debug2" , "cylindrical_debug8",]
fn = "ks_LaserEnergyProfile.dat"
l_hist = []

for r in run:
    filename = local_dir / r / fn
    if filename.exists():
        print(f"Reading {filename}")

    data = np.loadtxt(filename, comments="#")
    l = LaserEnergy(data)
    l_hist.append(l)

save_plot = True

# fig, (ax1, ax2) = plt.subplots(
#         2, 1,
#         figsize=(7, 6),
#         sharex=True
#     )
  
# ax1.plot(l_hist[0].time, l_hist[0].ein, lw=2, label="Energy in")
# ax1.plot(l_hist[0].time, l_hist[0].eout, lw=2, label="Energy out")
# ax1.plot(l_hist[0].time, l_hist[0].eabs, lw=2, label="Absorbed = in - out")

# ax1.set_ylabel("Energy [J]")
# ax1.set_title(run[0])
# ax1.grid(True, alpha=0.3)
# ax1.legend(frameon=False)

# ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

# ax2.plot(l_hist[1].time, l_hist[1].ein, lw=2)
# ax2.plot(l_hist[1].time, l_hist[1].eout, lw=2)
# ax2.plot(l_hist[1].time, l_hist[1].eabs, lw=2)

# ax2.set_xlabel("Time [ns]")
# ax2.set_ylabel("Energy [J]")
# ax2.set_title(run[1])
# ax2.grid(True, alpha=0.3)


fig, ax1 = plt.subplots(
        figsize=(7, 6),
        sharex=True
    )
  
ax1.plot(l_hist[0].time, l_hist[0].ein, lw=2, color = 'blue',label="Energy in, Cartesian")
ax1.plot(l_hist[0].time, l_hist[0].eout, lw=2, color = 'red', label="Energy out")
ax1.plot(l_hist[0].time, l_hist[0].eabs, lw=2, color = 'green', label="Absorbed = in - out")

ax1.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

ax1.plot(l_hist[1].time, l_hist[1].ein, lw=2, linestyle="--", color = 'blue', label="Cylindrical")
ax1.plot(l_hist[1].time, l_hist[1].eout, lw=2, linestyle="--", color = 'red')
ax1.plot(l_hist[1].time, l_hist[1].eabs, lw=2, linestyle="--", color = 'green')

ax1.set_xlabel("Time [ns]")
ax1.set_ylabel("Energy [J]")
ax1.set_title("Cartesian vs Cylindrical Laser Histories")
ax1.grid(True, alpha=0.3)
#ax2.set_ylim(-0.05, 1.05)
legend = ax1.legend(frameon=False, loc="upper left")
if save_plot:
    out = filename.with_name(filename.stem + "_energy_balance.png")
    plt.savefig(out, dpi=200, bbox_inches="tight")
    print(f"Saved {out}")

plt.show()
