import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from pathlib import Path
import yt
import glob
import os
import re
from flash_helpers import *


def plot_transverse_integrated_depo(
    ds,
    fp=None,
    geometry="cartesian",   # "cartesian" or "cylindrical"
    useMicrons=True,
    savePlots=False,
    saveDir=Path("."),
):
    """
    Plot transverse-integrated deposited energy.

    Cartesian:
        profile(x) = integral rho * depo dy

    Cylindrical R-Z:
        assumes x = r, y = z
        profile(z) = integral rho * depo * 2*pi*r dr

    Units:
        depo: erg/g
        rho:  g/cm^3

    Cartesian output:
        erg/cm^2

    Cylindrical output:
        erg/cm
    """

    # -----------------------
    # Hardcoded options
    # -----------------------
    xlim_um = [0, 1000]
    save_dpi = 150

    # -----------------------
    # Load data
    # -----------------------
    cg, dims = get_covering_grid(ds)
    sim_time_ns = get_sim_time_ns(ds)

    rho = np.asarray(get_field_array(cg, "flash", "dens")[:, :, 0])
    depo = np.asarray(get_field_array(cg, "flash", "depo")[:, :, 0])

    energy_density = rho * depo   # erg/cm^3

    nx, ny = energy_density.shape

    x0, x1 = float(ds.domain_left_edge[0]), float(ds.domain_right_edge[0])
    y0, y1 = float(ds.domain_left_edge[1]), float(ds.domain_right_edge[1])

    x_edges_cm = np.linspace(x0, x1, nx + 1)
    y_edges_cm = np.linspace(y0, y1, ny + 1)

    x_cent_cm = 0.5 * (x_edges_cm[:-1] + x_edges_cm[1:])
    y_cent_cm = 0.5 * (y_edges_cm[:-1] + y_edges_cm[1:])

    dx_cm = np.diff(x_edges_cm)
    dy_cm = np.diff(y_edges_cm)

    # -----------------------
    # Integrate transverse direction
    # -----------------------
    if geometry.lower() in ["cartesian", "cart"]:
        # x = transverse, y = longitudinal
        profile = np.nansum(
            energy_density * dx_cm[:, None],
            axis=0
        )

        coord_cm = y_cent_cm

        xlabel_cm = r"$y$ [cm]"
        xlabel_um = r"$y$ [$\mu$m]"
        ylabel = r"$\int \rho E_{\rm dep}\,dx$ [erg/cm$^2$]"
        title = rf"Cartesian transverse-integrated deposition, t = {sim_time_ns:.3f} ns"

    elif geometry.lower() in ["cylindrical", "cyl", "rz"]:
        # x = r, y = z
        r_cent_cm = x_cent_cm

        radial_weight = 2.0 * np.pi * r_cent_cm[:, None] * dx_cm[:, None]

        profile = np.nansum(energy_density * radial_weight, axis=0)

        coord_cm = y_cent_cm

        xlabel_cm = r"$z$ [cm]"
        xlabel_um = r"$z$ [$\mu$m]"
        ylabel = r"$\int \rho E_{\rm dep}\,2\pi r\,dr$ [erg/cm]"
        title = rf"Cylindrical radial-integrated deposition, t = {sim_time_ns:.3f} ns"
    else:
        raise ValueError("geometry must be 'cartesian' or 'cylindrical'.")

    if useMicrons:
        coord = coord_cm * 1e4
        xlabel = xlabel_um
    else:
        coord = coord_cm
        xlabel = xlabel_cm

    # -----------------------
    # Plot
    # -----------------------
    fig, ax = plt.subplots(figsize=(7, 4))

    ax.plot(coord, profile, lw=2)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.yaxis.get_offset_text().set_fontsize(11)

    if useMicrons and xlim_um is not None:
        ax.set_xlim(xlim_um)

    finite_profile = profile[np.isfinite(profile)]
    if finite_profile.size > 0:
        ymin = np.nanmin(finite_profile)
        ymax = np.nanmax(finite_profile)
        if ymin != ymax:
            pad = 0.05 * (ymax - ymin)
            ax.set_ylim(ymin - pad, ymax + pad)

    plt.tight_layout()

    if savePlots:
        saveDir = Path(saveDir)
        saveDir.mkdir(parents=True, exist_ok=True)

        if fp is None:
            fname = f"transverse_integrated_depo_{geometry}.png"
        else:
            fname = f"{Path(fp).stem}_transverse_integrated_depo_{geometry}.png"

        out = saveDir / fname
        plt.savefig(out, dpi=save_dpi, bbox_inches="tight")
        print(f"Saved {out}")

    return fig, ax, coord, profile



def plot_ray_power_transverse_profile(
    fp,
    sim_time_ns,
    geometry="cartesian",   # "cartesian" or "cylindrical"
    useMicrons=True,
    savePlots=False,
    saveDir=Path("."),
):
    """
    Plot laser ray power/intensity profile across transverse direction
    at a chosen axial location.

    ray_data columns assumed:
      0 = ray tag
      1 = x [cm]
      2 = y [cm]
      3 = z [cm]
      4 = power [erg/s]

    Cartesian:
      transverse = x
      axial      = y

    Cylindrical:
      transverse = x = r
      axial      = y = z-like simulation direction

    For cylindrical intensity:
      I(r) = P_bin / (2*pi*r*dr)

    For Cartesian:
      The code plots power per transverse length [W/cm],
      unless plot_intensity=False.
    """

    # -----------------------
    # Hardcoded options
    # -----------------------
    axial_pos_um = 100.0
    nbins = 32
    save_dpi = 150

    plot_intensity = True

    use_closest_point_per_ray = True

    # -----------------------
    # Load ray data
    # -----------------------
    ray_data = read_flash_rays(fp)

    tag_col = 0
    x_col = 1
    y_col = 2
    power_col = 4

    tags = ray_data[:, tag_col]

    if geometry.lower() in ["cartesian", "cart"]:
        transverse_label = r"$x$ [$\mu$m]" if useMicrons else r"$x$ [cm]"
        title_geom = "Cartesian"
        axial_label = "y"

    elif geometry.lower() in ["cylindrical", "cyl", "rz"]:
        transverse_label = r"$r$ [$\mu$m]" if useMicrons else r"$r$ [cm]"
        title_geom = "Cylindrical R-Z"
        axial_label = "z"

    else:
        raise ValueError("geometry must be 'cartesian' or 'cylindrical'.")

    axial_pos_cm = axial_pos_um * 1e-4

    # -----------------------
    # Pick ray power at chosen axial plane
    # -----------------------
    unique_tags = np.unique(tags)

    chosen_trans_cm = []
    chosen_power_W = []

    for tag in unique_tags:
        this_ray = ray_data[tags == tag]

        this_trans = this_ray[:, x_col]
        this_axial = this_ray[:, y_col]
        this_power_W = this_ray[:, power_col] * 1e-7

        if use_closest_point_per_ray:
            idx = np.argmin(np.abs(this_axial - axial_pos_cm))
            chosen_trans_cm.append(this_trans[idx])
            chosen_power_W.append(this_power_W[idx])

    chosen_trans_cm = np.asarray(chosen_trans_cm)
    chosen_power_W = np.asarray(chosen_power_W)

    # -----------------------
    # Coordinates for plotting
    # -----------------------
    if useMicrons:
        chosen_trans = chosen_trans_cm * 1e4
        axial_pos = axial_pos_um
        unit_label = r"$\mu$m"
    else:
        chosen_trans = chosen_trans_cm
        axial_pos = axial_pos_cm
        unit_label = "cm"

    # -----------------------
    # Bin ray power
    # -----------------------
    bin_power_W, bin_edges = np.histogram(
        chosen_trans,
        bins=nbins,
        weights=chosen_power_W
    )

    bin_counts, _ = np.histogram(
        chosen_trans,
        bins=bin_edges
    )

    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # -----------------------
    # Convert binned power to intensity-like quantity
    # -----------------------
    if plot_intensity:
        if useMicrons:
            bin_edges_cm = bin_edges * 1e-4
            bin_centers_cm = bin_centers * 1e-4
        else:
            bin_edges_cm = bin_edges
            bin_centers_cm = bin_centers

        dr_bins_cm = np.diff(bin_edges_cm)

        if geometry.lower() in ["cylindrical", "cyl", "rz"]:
            # Full axisymmetric annular area
            area_cm2 = 2.0 * np.pi * bin_centers_cm * dr_bins_cm

            yplot = np.full_like(bin_power_W, np.nan, dtype=float)
            good = area_cm2 > 0
            yplot[good] = bin_power_W[good] / area_cm2[good]

            y_label = r"Intensity [W/cm$^2$]"
            plot_label = "intensity"

        else:
            # Cartesian 2D has no actual out-of-plane beam area.
            # This is power per transverse length.
            dx_bins_cm = dr_bins_cm

            yplot = np.full_like(bin_power_W, np.nan, dtype=float)
            good = dx_bins_cm > 0
            yplot[good] = bin_power_W[good] / dx_bins_cm[good]

            y_label = r"Power per transverse length [W/cm]"
            plot_label = "power/length"

    else:
        yplot = bin_power_W
        y_label = "Summed ray power in bin [W]"
        plot_label = "power"

    # -----------------------
    # Plot
    # -----------------------
    fig, ax = plt.subplots(figsize=(7, 4))

    ax.step(bin_centers, yplot, where="mid", lw=2)

    ax.set_xlabel(transverse_label)
    ax.set_ylabel(y_label)
    ax.set_title(
        rf"{title_geom} ray {plot_label} profile, t = {sim_time_ns:.3f} ns"
        "\n"
        rf"{axial_label} = {axial_pos:.1f} {unit_label}"
    )

    ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax.grid(True, alpha=0.3)

    total_power_TW = np.nansum(chosen_power_W) / 1e12

    ax.text(
        0.98,
        0.95,
        f"rays = {len(unique_tags)}\n"
        f"sampled power = {total_power_TW:.3f} TW",
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )

    plt.tight_layout()

    if savePlots:
        saveDir = Path(saveDir)
        saveDir.mkdir(parents=True, exist_ok=True)

        safe_plot_label = plot_label.replace("/", "_per_").replace(" ", "_")

        out = saveDir / (
            f"{Path(fp).stem}_ray_{safe_plot_label}_profile_{geometry}.png"
        )

        plt.savefig(out, dpi=save_dpi, bbox_inches="tight")
        print(f"Saved {out}")
        
    return (
    fig,
    ax,
    bin_centers,
    yplot,
    bin_power_W,
    bin_counts,
    chosen_trans,
    chosen_power_W,
)
