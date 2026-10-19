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
from plotFLASH2d_profiles_cylindrical_only import *
from FLASH_ray_diagnostics import *
import os

eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi
c = constants.c
h = constants.h


# ============================================================
# User options
# ============================================================

# --------- Mac ------------
# runDir = baseDir / ".." / "Flash" / "2DCylindrical"/"Si3N4_cyl_5TW_1300psiGJS_cfl_0.1"
# file = "ks_hdf5_plt_cnt_0050"
# fp = runDir / file

# --------- Windows ------------
runDir = Path(r"C:\Simulation_data\FLASH\1D_spherical\iFP_bubble_T_7.5e10")  # directory containing the FLASH plot files
plt_by_time = False  # if True, plot by time using the flash_plotfile_times.txt file; if False, plot by file number

# Axis labels
useMicrons = True

# Save plots?
savePlots = True
saveDir = runDir / "plots"
saveDir.mkdir(parents=True, exist_ok=True)

# Ray option
rays = False
xlims = None
burnthrough_time_ns = None  # if None, will be determined automatically from the first file where deposition is nonzero at the chosen x location


# ============================================================
# Main loop over requested times
# ============================================================
i = 0
if plt_by_time:
    #target_time_ns = np.arange(0.6,0.8, 0.02)  # ns
    target_time_ns = [0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,2,3,4,5,6,7]  # ns
    shock_pos = np.zeros((2, len(target_time_ns)))
    shock_res = [] # shock position, shock time, shock jump ratio
    files_info = find_nearest_flash_file_from_table(target_time_ns, runDir, time_table="flash_plotfile_times.txt")
            
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

        #print_fields(ds, ad)   # comment if you don't want to see all the field names and units

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
        # field  = ["dens", "tele", "velx","pres"]
        # plot_1d_profiles_with_material_fractions(ds, fields = field, useMicrons=useMicrons, savePlots=savePlots, saveDir=saveDir, fp = fp, rays=ray_data)
        

        # fig, axes, burnthrough = plotFLASH1d_profiles(ds, field="dens", fp=fp, xlims=xlims, useMicrons=useMicrons, savePlots=savePlots, saveDir=saveDir, ray_data=ray_data)
        # # First file where deposition appears
        # if (burnthrough_time_ns is None and burnthrough["detected"]):

        #     burnthrough_time_ns = burnthrough["time_ns"]
        #     out = Path(saveDir) / "burnthrough_time.txt"
        #     with open(out, "w") as f:
        #         f.write(f"{burnthrough_time_ns:.10e}\n")

        #     print(f"Saved burn-through time to {out}")
        #     print("\n=================================")
        #     print(f"BURN THROUGH = " f"{burnthrough_time_ns:.6f} ns")
        #     print("=================================\n")
        #plt.close(fig)
        
        # ============================================================
        # Plot 2D
        # Comment/uncomment this block for 2D
        # ============================================================
        # field = "dens"  # field to plot in 2D
        # # # # ------- For the regular 2D plots --------------- #
        # # results = plot_2d_profiles_with_material_fractions(ds, fp, "flash", field, lineout_dir = "x", lineout_pos = 0, useMicrons = useMicrons, savePlot=savePlots, saveDir=saveDir,)

        # #results = plot_2d_ne_with_material_fractions(ds, fp, "flash", lineout_dir = "x", lineout_pos = 0, useMicrons = useMicrons, savePlot=savePlots, saveDir=saveDir,)

        # #results = plot_knudsen_number(ds, fp, cham={"A": 4, "Zion": 2}, targ={"A": 20.2, "Zion": 10}, gas={"A": 1.6, "Zion": 1.2}, test_species="gas", lineout_dir="x", lineout_pos=0.0, useMicrons=True, savePlot=True, saveDir=saveDir,cmap = "turbo", vmin=1e-7, vmax=1e1)

        # # -------- For the integrated energy deposition plots --------------- #  
        # # fig, ax, coord, profile = plot_transverse_integrated_depo(
        # # ds,
        # # fp=fp,
        # # geometry="cylindrical",
        # # savePlots=True,
        # # saveDir=saveDir,
        # # )

        # # -------- For the purple-green-red color scheme plots --------------- #  
        # fig, axes, shock_pos[0,i], shock_result = plotFLASH2d_profiles(ds, "flash", field, useMicrons, savePlots, saveDir, fp, rays=False)
        # shock_pos[1,i] = sim_time_ns
        # shock_res.append(shock_result)
        # i = i+1

        # ============================================================
        # Ray diagnostics
        # ============================================================

        if rays and ray_data is not None:
            ray_data = read_flash_rays(fp)
            fig, ax, bin_centers, yplot, bin_power_W, bin_counts, chosen_trans, chosen_power_W = plot_ray_power_transverse_profile(fp, sim_time_ns, geometry="cylindrical", savePlots=True, saveDir=saveDir)
            
            # plot number of rays if you want
            plt.figure()
            plt.step(bin_centers, bin_counts, where="mid")
            plt.xlabel("r [um]")
            plt.ylabel("Number of rays")
            plt.show()

            # diag = make_ray_diagnostic_figure(
            #     ds,
            #     fp,
            #     ray_data,
            #     geometry="cylindrical",
            #     savePlots=True,
            #     saveDir=saveDir,
            # )

else:
    fnums = np.arange(0,100,10)  # file numbers to plot
    files = find_flash_file_from_number(fnums, runDir)
    shock_res = []
    for i, file in enumerate(files):
        ds = yt.load(file)
        ad = ds.all_data()
        print_fields(ds, ad)   # comment if you don't want to see all the field names and units
        for field in ["tele","tion"]:
            *_, t, shock_result = plot_1d_field(ds, field=field, file_number=file, shock_jump_threshold = 5.0e-05, post_offset=0.01, savePlots=savePlots, saveDir=saveDir)
            # plot_2d_profiles_with_material_fractions(ds, file, "flash", field, lineout_dir = "x", lineout_pos = 1000, useMicrons = useMicrons, savePlot=savePlots, saveDir=saveDir,)

            if shock_result is not None:
                shock_res.append(shock_result) 


saveDir = Path(saveDir)
saveDir.mkdir(parents=True, exist_ok=True)

# # Save shock position vs time
# out = saveDir / "shock_position_vs_time.npy"
# np.save(out, shock_pos)
# print(f"Saved {out}")

# Save full shock-result dictionaries
out_res = saveDir / "shock_results.npy"
np.save(out_res, shock_res, allow_pickle=True)
print(f"Saved {out_res}")