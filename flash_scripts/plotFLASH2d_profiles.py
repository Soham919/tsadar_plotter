import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from pathlib import Path
import yt
from flash_helpers import *
from scipy import constants

eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi
c = constants.c
h = constants.h
A1 = 1 # H Gas
A2 = 28 # Si Targ
A3 = 4 # He cham

def plotFLASH2d_profiles(ds, ftype, field, useMicrons,
                  savePlots, saveDir, fp, rays=False):

    """
    Top: 2D FLASH field, optionally with ray overlay.
    Bottom: 1D lineout of ni, ne, and laser deposition.
    """

    # ============================================================
    # Hardcoded lineout settings
    # ============================================================
    lineout_axis = "y"   # "x" -> profile vs x at fixed y
                          # "y" -> profile vs y at fixed x

    lineout_pos_um = 0

    # ============================================================
    # Load data
    # ============================================================
    cg, dims = get_covering_grid(ds)
    sim_time_ns = get_sim_time_ns(ds)

    arr2d = get_field_array(cg, ftype, field)[:, :, 0]

    dens2d = get_field_array(cg, "flash", "dens")[:, :, 0]
    depo2d = get_field_array(cg, "flash", "depo")[:, :, 0]
    ye2d   = get_field_array(cg, "flash", "ye")[:, :, 0]

    gas2d  = get_field_array(cg, "flash", "gas")[:, :, 0]
    targ2d = get_field_array(cg, "flash", "targ")[:, :, 0]
    cham2d = get_field_array(cg, "flash", "cham")[:, :, 0]

    print(f"{field} 2D shape:", arr2d.shape)
    print(f"{field} min/max:", np.nanmin(arr2d), np.nanmax(arr2d))

    # ============================================================
    # Coordinates
    # ============================================================
    x0, x1 = float(ds.domain_left_edge[0]), float(ds.domain_right_edge[0])
    y0, y1 = float(ds.domain_left_edge[1]), float(ds.domain_right_edge[1])

    nx, ny = arr2d.shape

    x_edges_cm = np.linspace(x0, x1, nx + 1)
    y_edges_cm = np.linspace(y0, y1, ny + 1)

    x_centers_cm = 0.5 * (x_edges_cm[:-1] + x_edges_cm[1:])
    y_centers_cm = 0.5 * (y_edges_cm[:-1] + y_edges_cm[1:])

    if useMicrons:
        x_edges = x_edges_cm * 1e4
        y_edges = y_edges_cm * 1e4

        x_centers = x_centers_cm * 1e4
        y_centers = y_centers_cm * 1e4

        lineout_pos = lineout_pos_um

        xlabel_x = r"$x$ [$\mu$m]"
        xlabel_y = r"$y$ [$\mu$m]"
        coord_unit = r"$\mu$m"

    else:
        x_edges = x_edges_cm
        y_edges = y_edges_cm

        x_centers = x_centers_cm
        y_centers = y_centers_cm

        lineout_pos = lineout_pos_um * 1e-4

        xlabel_x = r"$x$ [cm]"
        xlabel_y = r"$y$ [cm]"
        coord_unit = "cm"

    X, Y = np.meshgrid(x_edges, y_edges, indexing="ij")

    # ============================================================
    # Extract 1D lineout
    # ============================================================
    if lineout_axis.lower() == "x":
        idx = np.argmin(np.abs(y_centers - lineout_pos))

        coord_plot = x_centers
        xlabel_line = xlabel_x

        dens = dens2d[:, idx]
        depo = depo2d[:, idx]
        ye   = ye2d[:, idx]

        gas  = gas2d[:, idx]
        targ = targ2d[:, idx]
        cham = cham2d[:, idx]

        line_location = y_centers[idx]
        line_title = rf"Lineout vs $x$ at $y = {line_location:.3f}$ {coord_unit}"

    elif lineout_axis.lower() == "y":
        idx = np.argmin(np.abs(x_centers - lineout_pos))

        coord_plot = y_centers
        xlabel_line = xlabel_y

        dens = dens2d[idx, :]
        depo = depo2d[idx, :]
        ye   = ye2d[idx, :]

        gas  = gas2d[idx, :]
        targ = targ2d[idx, :]
        cham = cham2d[idx, :]

        line_location = x_centers[idx]
        line_title = rf"Lineout vs $y$ at $x = {line_location:.3f}$ {coord_unit}"

    else:
        raise ValueError("lineout_axis must be either 'x' or 'y'")

    print(f"Lineout axis: {lineout_axis}")
    print(f"Requested lineout position: {lineout_pos_um} microns")
    print(f"Nearest lineout location used: {line_location:.4f} {coord_unit}")

    # ============================================================
    # Convert lineout fields
    # ============================================================
    mp_g = mp * 1000

    gas_dens    = dens * gas  / (A1 * mp_g)
    target_dens = dens * targ / (A2 * mp_g)
    cham_dens   = dens * cham / (A3 * mp_g)

    ni = gas_dens + target_dens + cham_dens
    ne = ye * dens / mp_g

    las_depo = dens * depo * 1e-7

    # ============================================================
    # Rays for 2D overlay only
    # ============================================================
    if rays:
        ray_data = read_flash_rays(fp)
    else:
        ray_data = None

    # ============================================================
    # Make figure
    # ============================================================
    fig, (ax2d, ax1d) = plt.subplots(2,1,figsize=(9, 7),
    gridspec_kw={"height_ratios": [2.2, 1]})

    # ============================================================
    # Top: 2D plot
    # ============================================================
    positive = arr2d[arr2d > 0]

    if positive.size > 0:
        norm = LogNorm(vmin=np.nanmin(positive), vmax=np.nanmax(arr2d))
    else:
        print(f"Warning: {field} has no positive values. Using linear scale.")
        norm = Normalize(vmin=np.nanmin(arr2d), vmax=np.nanmax(arr2d))

    # Rotate top 2D plot 90 degrees clockwise
    # old x -> new y
    # old y -> new x, flipped
    arr2d_rot = arr2d.T[:, ::-1]
    Xrot, Yrot = np.meshgrid(y_edges, x_edges[::-1],indexing="ij")
    pcm = ax2d.pcolormesh(Xrot, Yrot, arr2d_rot, shading="auto", cmap="plasma", norm=norm)

    if rays and ray_data is not None:
        plot_rays(
            ax2d,
            ray_data,
            xscale=1e4,
            yscale=1e4,
            color="w",
            lw=0.8
        )

    if lineout_axis.lower() == "x":
    # original horizontal y-line becomes vertical line after clockwise rotation
        ax2d.axvline(line_location, color="w", linestyle="--", linewidth=1.2)
    else:
    # original vertical x-line becomes horizontal line after clockwise rotation
        ax2d.axhline(line_location, color="w", linestyle="--", linewidth=1.2)

    ax2d.set_xlabel(xlabel_y)
    ax2d.set_ylabel(xlabel_x)
    ax2d.set_title(f"{field}, t = {sim_time_ns:.3f} ns")
    ax2d.set_aspect("equal")

    cbar = fig.colorbar(pcm, ax=ax2d, orientation="horizontal", location="top", fraction=0.08, pad=0.12, aspect=30)
    cbar.set_label(field)

    # ============================================================
    # Bottom: 1D hydro + laser deposition lineout only
    # ============================================================
    c0, c1 = ["green", "red"]

    l1, = ax1d.plot(coord_plot, ni, lw=2, color=c0, label=r"$n_i$")
    l2, = ax1d.plot(coord_plot, ne, lw=2, color=c1, label=r"$n_e$")

    ax1d.set_ylabel(r"$n$ [$cm^{-3}$]")
    ax1d.set_xlabel(xlabel_line)
    ax1d.set_title(line_title)
    ax1d.grid(True, alpha=0.3)

    ax1d_2 = ax1d.twinx()

    c2 = "blueviolet"

    l3, = ax1d_2.plot(
        coord_plot,
        las_depo,
        color=c2,
        linewidth=2,
        linestyle="--",
        label=r"$E_{dep}$"
    )

    ax1d_2.set_ylabel(r"$E_{dep}$ [J/cm$^3$]", color=c2)
    ax1d_2.tick_params(axis="y", colors=c2)
    ax1d_2.spines["right"].set_color(c2)
    ax1d_2.spines["right"].set_linewidth(1.5)

    lines = [l1, l2, l3]
    labels = [line.get_label() for line in lines]
    ax1d.legend(lines, labels, loc="upper right", frameon=False)

    plt.tight_layout()

    if savePlots:
        out = saveDir / f"{Path(fp).stem}_{field}_2D_with_lineout.png"
        plt.savefig(out, dpi=200)
        print(f"Saved {out}")

    plt.show()

    return fig, (ax2d, ax1d, ax1d_2)