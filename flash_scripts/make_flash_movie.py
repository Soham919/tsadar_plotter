import subprocess
from attrs import field
import numpy as np
import matplotlib.pyplot as plt
import yt
from flash_helpers import *
from flash_helpers import _auto_downsample_2d
from scipy import constants
from pathlib import Path
from matplotlib.colors import LogNorm, Normalize

eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi
c = constants.c
h = constants.h

# ============================================================
# GLOBAL PLOT STYLE
# ============================================================

# ---- Fonts ----
FONT_FAMILY = "serif"
FONT_SIZE = 24
LABEL_SIZE = 26
TITLE_SIZE = 28
TICK_SIZE = 24
LEGEND_SIZE = 22

# ---- Lines / axes ----
AXIS_WIDTH = 1.2
TICK_WIDTH = 1.2
TICK_LENGTH = 5
LINE_WIDTH = 2.0

plt.rcParams.update({
    # Font
    "font.family": FONT_FAMILY,
    "font.size": FONT_SIZE,

    # Use Computer Modern-style math
    "mathtext.fontset": "cm",

    # Axes
    "axes.labelsize": LABEL_SIZE,
    "axes.titlesize": TITLE_SIZE,
    "axes.linewidth": AXIS_WIDTH,

    # Ticks
    "xtick.labelsize": TICK_SIZE,
    "ytick.labelsize": TICK_SIZE,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.major.width": TICK_WIDTH,
    "ytick.major.width": TICK_WIDTH,
    "xtick.major.size": TICK_LENGTH,
    "ytick.major.size": TICK_LENGTH,

    # Legend
    "legend.fontsize": LEGEND_SIZE,

    # Plot lines
    "lines.linewidth": LINE_WIDTH,
})



def make_flash_movie(files, ftype, field, output="flash_movie.mp4", lineout_pos=0.0, fps=10,
                     useMicrons=True, saveDir=Path("."), max_plot_pixels=1200,
                     vmin=None, vmax=None, cmap="viridis"):
    """
    Make a FLASH movie.

    Top:    2D map of the selected field.
    Bottom: lineout of the same field along y at fixed x=lineout_pos,
            colored by dominant material.

    Shock detection, if enabled below, always uses the density lineout.
    """

    # ==========================================================
    # SETTINGS TO CHANGE HERE
    # ==========================================================
    SHOCK_DETECTION = True
    jump_threshold, shock_window, shock_mode = 1e-5, 5, "up"
    shock_min_coord, shock_max_coord = None, None

    movie_xlim, movie_ylim = None, None       # Example: (4500,6500), (900,1100)
    lineout_ylim = None                       # Example: (0,1000)
    lineout_logscale = False
    material_threshold = 0.5

    FORCE_LINEAR = field.lower() in ["velx", "vely", "velz"]
    MAKE_FRAMES = True  # if True, make frames and save as PNG; if False, skip frame generation and just make movie from existing frames

    # ==========================================================
    # Setup
    # ==========================================================
    files, saveDir = [Path(f) for f in files], Path(saveDir)
    if not files: raise ValueError("No FLASH plotfiles supplied.")

    saveDir.mkdir(parents=True, exist_ok=True)
    frameDir = saveDir / f"{Path(output).stem}_frames"
    frameDir.mkdir(parents=True, exist_ok=True)

    #coordinate_scale = 1e-4 if useMicrons else 1.0
    #coordinate_unit = r"$\mu$m" if useMicrons else "cm"

    coordinate_scale = 10 if useMicrons else 1.0
    coordinate_unit = r"mm" if useMicrons else "cm"
    

    # ==========================================================
    # Load one frame
    # ==========================================================
    def load_frame(fp):
        ds = yt.load(str(fp))
        cg, dims = get_covering_grid(ds)
        time = get_sim_time_ns(ds)

        arr = np.asarray(get_field_array(cg, ftype, field)[:, :, 0])
        dens = np.asarray(get_field_array(cg, ftype, "dens")[:, :, 0])
        cham = np.asarray(get_field_array(cg, ftype, "cham")[:, :, 0])
        gas = np.asarray(get_field_array(cg, ftype, "gas")[:, :, 0])
        targ = np.asarray(get_field_array(cg, ftype, "targ")[:, :, 0])
        tele = np.asarray(get_field_array(cg, ftype, "tele")[:, :, 0])
        tion = np.asarray(get_field_array(cg, ftype, "tion")[:, :, 0])
        vely = np.asarray(get_field_array(cg, ftype, "vely")[:, :, 0])

        if field.lower() in ["tion", "tele"]: arr = arr * kb/e
        if field.lower() in ["dens"]: arr *= 1000

        nx, ny = arr.shape
        xe = np.linspace(float(ds.domain_left_edge[0]), float(ds.domain_right_edge[0]), nx+1) * coordinate_scale
        ye = np.linspace(float(ds.domain_left_edge[1]), float(ds.domain_right_edge[1]), ny+1) * coordinate_scale
        xc, yc = 0.5*(xe[:-1]+xe[1:]), 0.5*(ye[:-1]+ye[1:])

        ix = int(np.argmin(np.abs(xc-lineout_pos)))

        return {
            "time": time, "field2d": arr, "field_line": arr[ix, :], "dens_line": dens[ix, :], "tele_line": tele[ix, :], "tion_line": tion[ix, :], "vely_line": vely[ix, :],
            "cham": cham[ix, :], "gas": gas[ix, :], "targ": targ[ix, :],
            "xe": xe, "ye": ye, "yc": yc,  "line_location": xc[ix], "sim_time" : time,
        }

    # ==========================================================
    # Find one global color scale for entire movie
    # ==========================================================
    global_vmin, global_vmax = vmin, vmax

    if MAKE_FRAMES and (global_vmin is None or global_vmax is None):
        print("Finding global color limits...")
        all_min, all_max = np.inf, -np.inf

        for i, fp in enumerate(files):
            print(f"Scanning {i+1}/{len(files)}: {fp.name}")
            arr = load_frame(fp)["field2d"]
            finite = arr[np.isfinite(arr)]
            if not FORCE_LINEAR: finite = finite[finite > 0]
            if finite.size == 0: continue
            all_min, all_max = min(all_min, np.min(finite)), max(all_max, np.max(finite))

        if global_vmin is None: global_vmin = all_min
        if global_vmax is None: global_vmax = all_max

    if MAKE_FRAMES:
        print(f"{field}: {global_vmin:.4g} → {global_vmax:.4g}")
        norm = Normalize(vmin=global_vmin, vmax=global_vmax) if FORCE_LINEAR else LogNorm(vmin=global_vmin, vmax=global_vmax)

    # ==========================================================
    # Make frames
    # ==========================================================
    if MAKE_FRAMES:
        for n, fp in enumerate(files):
            print(f"Frame {n+1}/{len(files)}: {fp.name}")
            d = load_frame(fp)

            arr, line, dens, tele, tion, vely = d["field2d"], d["field_line"], d["dens_line"], d["tele_line"], d["tion_line"], d["vely_line"]
            cham, gas, targ, yc = d["cham"], d["gas"], d["targ"], d["yc"]
            sim_time_ns = d["sim_time"]
            image, stride = _auto_downsample_2d(arr, max_pixels=max_plot_pixels)
            extent = [d["ye"][0], d["ye"][-1], d["xe"][0], d["xe"][-1]]

            # Dominant material: 0=target, 1=gas, 2=chamber
            fractions = np.vstack((targ, gas, cham))
            material = np.argmax(fractions, axis=0)
            strength = np.max(fractions, axis=0)

            target = (material == 0) & (strength >= material_threshold)
            gas_mask = (material == 1) & (strength >= material_threshold)
            chamber = (material == 2) & (strength >= material_threshold)
            mixed = strength < material_threshold

            # Shock detection still uses density
            shock = None
            if SHOCK_DETECTION:

                shock_result = detect_first_density_jump_right_to_left(
                                yc, sim_time_ns, dens, tele, tion, vely,
                                gas, targ, cham, jump_threshold = jump_threshold,
                                window=shock_window,
                                post_width=500,
                                post_offset=50,
                                min_coord=shock_min_coord,
                                max_coord=shock_max_coord,
                                statistic="median",
                                )
            # ======================================================
            # Figure
            # ======================================================

            TARGET_COLOR = "tab:blue"
            GAS_COLOR = "tab:orange"
            CHAMBER_COLOR = "tab:green"
            MIXED_COLOR = "orange"

            fig, (ax2d, axline) = plt.subplots(2, 1, figsize=(6.5, 8), sharex=True, gridspec_kw={"height_ratios": [3, 1.25],"hspace": 0.08})

            # Fixed margins -> identical geometry for every frame
            fig.subplots_adjust(left=0.18, right=0.95, bottom=0.11, top=0.88, hspace=0.08)
            im = ax2d.imshow(image, origin="lower", extent=extent, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm, rasterized=True)
            ax2d.set_ylabel(rf"$x$ ({coordinate_unit})")
            #ax2d.set_title(rf"$\rho$, $t={d['time']:.3f}$ ns")

            ax2d.axhline(d["line_location"], color="white", linestyle="--", linewidth=2.0, alpha=0.8)

            # ======================================================
            # Horizontal colorbar above the 2D plot
            # ======================================================
            cbar = fig.colorbar(im, ax=ax2d, orientation="horizontal", location="top", fraction=0.06, pad=0.04, aspect=30)
            cbar.ax.tick_params(width=TICK_WIDTH,length=TICK_LENGTH, labelsize=TICK_SIZE, pad=3)
            cbar.set_label(r"$\rho$ (kg/m$^3$)", fontsize=LABEL_SIZE, labelpad=10)
            cbar.outline.set_linewidth(AXIS_WIDTH)

            # ======================================================
            # Same-field lineout, colored by material
            # ======================================================

            axline.plot(yc, np.where(target, line, np.nan), label="Target", color=TARGET_COLOR)
            axline.plot(yc, np.where(gas_mask, line, np.nan), label="Gas", color=GAS_COLOR)
            axline.plot(yc, np.where(chamber, line, np.nan), label="Chamber", color=CHAMBER_COLOR)
            if np.any(mixed): axline.plot(yc, np.where(mixed, line, np.nan), linestyle="--", label="Mixed", color=MIXED_COLOR)

            axline.set_xlabel(rf"$y$ ({coordinate_unit})")
            axline.set_ylabel(r"$\rho$ (kg/m$^3$)")
            axline.grid(True, alpha=0.25)
            #axline.legend(loc="best", frameon=False, ncol=4)

            if lineout_logscale: axline.set_yscale("log")
            if lineout_ylim is not None: axline.set_ylim(*lineout_ylim)

            # ======================================================
            # Shock marker
            # ======================================================
            if SHOCK_DETECTION and shock_result["found"]:
                shock_pos = shock_result["shock_coord"]
                ax2d.axvline(shock_pos, color="white", linestyle="--", linewidth=2)
                axline.axvline(shock_pos, color="black", linestyle="--", linewidth=1.8)
            # ======================================================
            # Plot limits
            # ======================================================
            xlimits = movie_xlim if movie_xlim is not None else (d["ye"][0], d["ye"][-1])
            ax2d.set_xlim(*xlimits)

            if movie_ylim is not None: ax2d.set_ylim(*movie_ylim)

            # ======================================================
            # Save frame
            # ======================================================
            frame_path = frameDir / f"frame_{n:05d}.pdf"    # png for movie, pdf for plots
            fig.savefig(frame_path)
            plt.close(fig)

    # ==========================================================
    # Turn PNG frames into MP4
    # ==========================================================
    output_path = saveDir / output
    FFMPEG = r"C:\Users\sban\AppData\Local\Microsoft\WinGet\Packages\Gyan.FFmpeg_Microsoft.Winget.Source_8wekyb3d8bbwe\ffmpeg-9.0.2-full_build\bin\ffmpeg.exe"

    command = [FFMPEG, "-y", "-framerate", str(fps),
           "-i", str(frameDir / "frame_%05d.png"),
           "-vf", "pad=ceil(iw/2)*2:ceil(ih/2)*2",
           "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18",
           str(output_path)]


    try:
        subprocess.run(command, check=True)
    except FileNotFoundError:
        raise RuntimeError("ffmpeg was not found. Install ffmpeg and make sure it is available from your terminal.")

    print(f"Movie saved to: {output_path}")
    return output_path

# ============================================================
# User options
# ============================================================

# --------- Mac ------------
# runDir = baseDir / ".." / "Flash" / "2DCylindrical"/"Si3N4_cyl_5TW_1300psiGJS_cfl_0.1"
# file = "ks_hdf5_plt_cnt_0050"
# fp = runDir / file

# --------- Windows ------------
runDir = Path(r"C:\Simulation_data\FLASH\2D\FV4.7.1_carbuncle_debugs\26A_GJS_5TW")  # directory containing the FLASH plot files
plt_by_time = True  # if True, plot by time using the flash_plotfile_times.txt file; if False, plot by file number

# Axis labels
useMicrons = True

# Save plots?
savePlots = True
saveDir = runDir / "plots"
saveDir.mkdir(parents=True, exist_ok=True)

# ============================================================
# Main loop over requested times
# ============================================================
i = 0
if plt_by_time:
    #target_time_ns = np.arange(0.6,0.8, 0.02)  # ns
    target_time_ns = np.arange(0.0, 10, 1) # ns
    shock_pos = np.zeros((2, len(target_time_ns)))
    files_info = find_nearest_flash_file_from_table(target_time_ns, runDir, time_table="flash_plotfile_times.txt")

field = "dens"  # field to plot

make_flash_movie(field=field, files=[f[0] for f in files_info], ftype="flash", output=f"inferno_{field}.mp4",
                 lineout_pos=0.0, fps=8, useMicrons=useMicrons, saveDir=saveDir, max_plot_pixels=1200,
                 vmin=1e-3, vmax=1, cmap="inferno" )
