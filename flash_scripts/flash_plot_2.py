import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from pathlib import Path
import yt
from flash_helpers import *
from scipy import constants
from plotFLASH1d_profiles import *
from plotFLASH2d_profiles import *
import os

eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi
c = constants.c
h = constants.h
A1 = 1 # H
A2 = 28 # Si

# ============================================================
# User options
# ============================================================

baseDir = Path().resolve().parent

# --------- Mac ------------
runDir = baseDir / ".." / "Flash" / "test_runs" / "10um_spot_2"
# file = "ks_hdf5_plt_cnt_0050"
# fp = runDir / file

# --------- Windows ------------
runDir = baseDir / "FLASH" / "1D" / "1mm_spot" / "1D_1um_resolve"

target_time_ns = [0.1, 0.2, 0.3, 0.31, 0.32, 0.33, 0.34, 0.35, 0.4, 0.5, 0.6, 0.7, 0.8, 1,
1.5, 2,2.5, 3,3.5, 4,4.5, 5,5.5, 6,6.5, 7,7.5, 8,8.5, 9]  # ns

files_info = find_nearest_flash_file_from_table(
    target_time_ns,
    runDir,
    time_table="flash_plotfile_times.txt"
)

# Axis labels
useMicrons = True

# Save plots?
savePlots = True
saveDir = runDir / "plots"
saveDir.mkdir(parents=True, exist_ok=True)

# Ray option
rays = False
xlims = None


# ============================================================
# Main loop over requested times
# ============================================================

# ============================================================
# Main loop over requested times
# ============================================================

for filename, file_num, matched_time_s, matched_time_ns in files_info:

    fp = filename

    print("\n" + "=" * 70)
    print(f"Loading file number {file_num}")
    print(f"File: {fp}")
    print(f"Matched time = {matched_time_ns:.4f} ns")
    print("=" * 70)

    ds = yt.load(fp)
    ad = ds.all_data()

    sim_time_ns = get_sim_time_ns(ds)

    print("current_time =", ds.current_time)
    print(f"time = {sim_time_ns:.4f} ns")

    print_fields(ds, ad)

    cg, dims = get_covering_grid(ds)

    # ============================================================
    # Read rays if needed
    # ============================================================
    if rays:
        ray_data = read_flash_rays(fp)
    else:
        ray_data = None

    # ============================================================
    # Plot 1D
    # Comment/uncomment this block for 1D
    # ============================================================

    # fig, axes = plotFLASH1d_profiles(
    #     ds,
    #     cg,
    #     dims,
    #     xlims=xlims,
    #     ray_data=ray_data
    # )

    # if savePlots:
    #     save_name = saveDir / f"FLASH_1D_profiles_t_{sim_time_ns:.3f}_ns.png"
    #     fig.savefig(save_name, dpi=300, bbox_inches="tight")
    #     print(f"Saved: {save_name}")

    # plt.close(fig)

    # ============================================================
    # Plot 2D
    # Comment/uncomment this block for 2D
    # ============================================================

    # plot_2d_field(
    #     ds,
    #     "flash",
    #     field2D,
    #     useMicrons,
    #     savePlots,
    #     saveDir,
    #     fp,
    #     rays
    # )
    field = "dens"
    plotFLASH2d_profiles(ds, "flash", field, useMicrons, savePlots, saveDir, fp, rays)


    # ============================================================
    # Ray diagnostics
    # ============================================================

    if rays and ray_data is not None:
        tags = ray_data[:, 0]

        unique_tags, counts = np.unique(tags, return_counts=True)

        print("Number of rays:", len(unique_tags))

        ray_initial_power_W = []

        for tag in unique_tags:
            r = ray_data[tags == tag]
            ray_initial_power_W.append(r[0, 4] * 1e-7)

        ray_initial_power_W = np.array(ray_initial_power_W)

        print("Total initial laser power TW:", ray_initial_power_W.sum() / 1e12)
        print(
            "Initial power per ray W min/max:",
            ray_initial_power_W.min(),
            ray_initial_power_W.max()
        )
# ============================================================
# Plot everything one by one
# ============================================================
# Plot mode:
#   "auto" -> detect 1D or 2D
#   "1d"   -> force 1D profile plots
#   "2d"   -> force 2D pcolormesh

#   plotMode = "auto"
        
# fields1D = [
#     "dens",
#     "depo",
#     "tele",
#     "tion",
#     "trad",
#     "cham",
#     "targ",
#     "gas",
# ]
        
# # For 2D plots
# field2D = "dens"

# type = "gas"
# if plotMode.lower() == "auto":
#     detectedMode = detect_plot_dim_from_ds(ds)
# else:
#     detectedMode = plotMode.lower()

# print("plot mode:", detectedMode)

# if detectedMode == "1d":
#     plot_1d_profiles(ds, cg, fields1D, sim_time_ns,useMicrons, savePlots, saveDir)

# elif detectedMode == "2d":
#     plot_2d_field(ds, cg, type, field2D, sim_time_ns,useMicrons, savePlots, saveDir)

# else:
#     raise ValueError(f"Unsupported detected plot mode: {detectedMode}")