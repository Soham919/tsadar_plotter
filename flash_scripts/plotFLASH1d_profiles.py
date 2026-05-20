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

def plotFLASH1d_profiles(ds, cg, dims=None, xlims=None, useMicrons=True, ray_data=None):
    sim_time_ns = get_sim_time_ns(ds)

    # ---- get raw data ---- #
    dens = squeeze_to_1d(get_field_array(cg, "flash", "dens"))
    pres = squeeze_to_1d(get_field_array(cg, "flash", "pres"))
    gas  = squeeze_to_1d(get_field_array(cg, "flash", "gas"))
    targ = squeeze_to_1d(get_field_array(cg, "flash", "targ"))
    cham = squeeze_to_1d(get_field_array(cg, "flash", "cham"))
    depo = squeeze_to_1d(get_field_array(cg, "flash", "depo"))
    ye   = squeeze_to_1d(get_field_array(cg, "flash", "ye"))

    # ---- convert ---- #
    mp_g = mp * 1000  # kg -> g, if mp is in kg
    gas_dens    = dens * gas  / (A1 * mp_g)
    target_dens = dens * targ / (A2 * mp_g)
    cham_dens   = dens * cham / (A3 * mp_g)

    ni = gas_dens + target_dens + cham_dens
    ne = ye * dens / mp_g

    las_depo = dens * depo * 1e-7  # erg/cm^3 -> J/cm^3

    # ---- x axis ---- #
    nx = gas.size
    x0 = float(ds.domain_left_edge[0])
    x1 = float(ds.domain_right_edge[0])

    x_cm = np.linspace(x0, x1, nx, endpoint=False)
    dx = (x1 - x0) / nx
    x_cm = x_cm + 0.5 * dx
    x_plot, xlabel = coord_to_plot_units(x_cm, useMicrons)

    # ---- quick check: is laser deposition nonzero near chosen x? ---- #
    check_x_um = 3600.0   # choose x location in microns
    tol = 1e-30          # tolerance for "nonzero"

    RED = "\033[91m"
    RESET = "\033[0m"

    if useMicrons:
        x_for_check = x_plot
    else:
        x_for_check, _ = coord_to_plot_units(x_cm, True)

    ix_check = np.argmin(np.abs(x_for_check - check_x_um))
    val_check = las_depo[ix_check]

    if abs(val_check) > tol:
        print(RED +
            f"Laser deposition is NONZERO near x = {check_x_um:.3g} um "
            f"(nearest x = {x_for_check[ix_check]:.3g} um, value = {val_check:.3e} J/cm^3)"
            + RESET)
    else:
        print(RED +
            f"Laser deposition is zero near x = {check_x_um:.3g} um "
            f"(nearest x = {x_for_check[ix_check]:.3g} um, value = {val_check:.3e} J/cm^3)"
            + RESET)
    #------------------------------------------------------------------ #

    # Decide plot x-limits once, using hydro x-axis if xlims is not passed
    if xlims is None:
        plot_xlims = (np.nanmin(x_plot), np.nanmax(x_plot))
    else:
        plot_xlims = xlims

    # ---- figure ---- #
    fig, (ax, ax_ray) = plt.subplots(
        2, 1,
        figsize=(9, 6),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1]}
    )

    # =========================
    # Top subplot: density + Edep
    # =========================
    c0, c1 = ["green", "red"]

    l1, = ax.plot(x_plot, ni, lw=2, color=c0, label=r"$n_i$")
    l2, = ax.plot(x_plot, ne, lw=2, color=c1, label=r"$n_e$")

    ax.set_ylabel(r"$n$ [$cm^{-3}$]")
    ax.set_xlabel(xlabel)
    ax.set_title(f"FLASH 1D Density + Laser, t = {sim_time_ns:.3f} ns")
    ax.grid(True, alpha=0.3)
    ax.set_xlim(plot_xlims)

    ax2 = ax.twinx()
    c2 = "blueviolet"

    l3, = ax2.plot(
        x_plot, las_depo,
        color=c2,
        linewidth=2,
        linestyle="--",
        label=r"$E_{dep}$"
    )

    ax2.set_ylabel(r"$E_{dep}$ [J/cm$^3$]", color=c2)
    #ax2.set_ylim([0,10])
    ax2.tick_params(axis="y", colors=c2)
    ax2.spines["right"].set_color(c2)
    ax2.spines["right"].set_linewidth(1.5)
    ax2.set_xlim(plot_xlims)

    lines = [l1, l2, l3]
    labels = [line.get_label() for line in lines]
    ax.legend(lines, labels, loc="upper right", frameon=False)

    # =========================
    # Bottom subplot: ray power
    # =========================

    # Treat empty arrays the same as no ray_data
    has_ray_data = ray_data is not None and len(ray_data) > 0

    if has_ray_data:
        tags = ray_data[:, 0].astype(int)
        unique_tags = np.unique(tags)

        for tag in unique_tags:
            r = ray_data[tags == tag]

            ray_x_cm = r[:, 1]            # column 1 = x
            ray_power_W = r[:, 4] * 1e-7  # erg/s -> W

            # sort by x-position
            order = np.argsort(ray_x_cm)
            ray_x_cm = ray_x_cm[order]
            ray_power_W = ray_power_W[order]

            ray_x_plot, _ = coord_to_plot_units(ray_x_cm, useMicrons)

            ax_ray.plot(
                ray_x_plot,
                ray_power_W,
                lw=2,
                marker="o",
                ms=3,
                label=f"ray {tag}"
            )

        if len(unique_tags) > 1:
            ax_ray.legend(frameon=False, loc="best")

    else:
        ax_ray.text(
            0.5, 0.5,
            "No Ray Data passed",
            transform=ax_ray.transAxes,
            ha="center",
            va="center"
        )

        # Optional: keep bottom panel from looking totally empty
        ax_ray.set_ylim(0, 1)

    ax_ray.set_ylabel("Ray power [W]")
    ax_ray.set_xlabel(xlabel)
    ax_ray.grid(True, alpha=0.3)

    # Important: force shared x-axis limits even when there is no ray data
    ax_ray.set_xlim(plot_xlims)

    # Make sure x tick labels show on the top plot too, despite sharex=True
    ax.tick_params(axis="x", labelbottom=True)

    # ax.set_ylim([0, 1e21])
    # ax2.set_ylim([0, 1])
    fig.tight_layout()
    plt.show()

    return fig, (ax, ax2, ax_ray)
