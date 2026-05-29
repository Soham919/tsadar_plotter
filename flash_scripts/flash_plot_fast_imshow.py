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


def _auto_downsample_2d(arr, max_pixels=1200):
    """
    Downsample a 2D array for display only.

    The lineout still uses the full-resolution arrays.
    max_pixels controls the maximum displayed size along the largest axis.
    """
    n0, n1 = arr.shape
    stride = max(1, int(np.ceil(max(n0, n1) / max_pixels)))
    return arr[::stride, ::stride], stride


def plotFLASH2d_profiles(ds, ftype, field, useMicrons,
                  savePlots, saveDir, fp, rays=False):

    """
    Top: 2D FLASH field, optionally with ray overlay.
    Bottom: 1D lineout of ni, ne, and laser deposition.

    Speed changes:
      - Uses imshow instead of pcolormesh for the 2D image.
      - Avoids creating large 2D coordinate meshgrids.
      - Downsamples only the displayed 2D image if it is very large.
      - Does not load tele/tion unless the temperature-lineout toggle is enabled.
    """

    # ============================================================
    # Hardcoded settings
    # ============================================================
    lineout_axis = "y"   # "x" -> profile vs x at fixed y
                          # "y" -> profile vs y at fixed x

    lineout_pos_um = 0

    # Display speed knob:
    # If your image is huge, only the displayed 2D image is downsampled.
    # The 1D lineout is still extracted from full-resolution data.
    max_plot_pixels = 1200

    # imshow is much faster than pcolormesh for regularly spaced grids.
    use_fast_imshow = True

    # Keep this False unless you actually want the Te/Ti lineout section.
    # Loading tele/tion can be expensive for large FLASH files.
    plot_temperature_lineout = False

    # Save speed/size knob
    save_dpi = 150

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

    if plot_temperature_lineout:
        tele2d = get_field_array(cg, "flash", "tele")[:, :, 0]
        tion2d = get_field_array(cg, "flash", "tion")[:, :, 0]
    else:
        tele2d = None
        tion2d = None

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

        ray_scale = 1e4

    else:
        x_edges = x_edges_cm
        y_edges = y_edges_cm

        x_centers = x_centers_cm
        y_centers = y_centers_cm

        lineout_pos = lineout_pos_um * 1e-4

        xlabel_x = r"$x$ [cm]"
        xlabel_y = r"$y$ [cm]"
        coord_unit = "cm"

        ray_scale = 1.0

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

        if plot_temperature_lineout:
            tele = tele2d[:, idx]
            tion = tion2d[:, idx]

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

        if plot_temperature_lineout:
            tele = tele2d[idx, :]
            tion = tion2d[idx, :]

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
    fig, (ax2d, ax1d) = plt.subplots(
        2, 1,
        figsize=(9, 7),
        gridspec_kw={"height_ratios": [2.2, 1]},
        constrained_layout=False
    )

    # ============================================================
    # Top: 2D plot
    # ============================================================
    positive = arr2d[arr2d > 0]

    # clims if you want to set manually:
    # vmin = 5e-4
    # vmax = 5e-2
    # norm = LogNorm(vmin=vmin, vmax=vmax)

    if positive.size > 0:
        norm = LogNorm(vmin=np.nanmin(positive), vmax=np.nanmax(arr2d))
    else:
        print(f"Warning: {field} has no positive values. Using linear scale.")
        norm = Normalize(vmin=np.nanmin(arr2d), vmax=np.nanmax(arr2d))

    # Rotate top 2D plot 90 degrees clockwise:
    # horizontal axis = old y
    # vertical axis   = old x
    #
    # For imshow, array shape is [vertical, horizontal], so use:
    # rows    = old x, flipped so old x increases upward
    # columns = old y
    arr_img = arr2d[::-1, :]

    if use_fast_imshow:
        arr_img_plot, display_stride = _auto_downsample_2d(arr_img, max_pixels=max_plot_pixels)

        if display_stride > 1:
            print(
                f"Display downsampling enabled: showing every {display_stride}th cell "
                f"for 2D image only. Lineout remains full resolution."
            )

        pcm = ax2d.imshow(
            arr_img_plot,
            extent=[y_edges[0], y_edges[-1], x_edges[0], x_edges[-1]],
            origin="upper",
            aspect="equal",
            interpolation="nearest",
            cmap="plasma",
            norm=norm
        )

    else:
        # Slower fallback. Useful only if you truly need pcolormesh cell edges.
        arr2d_rot = arr2d.T[:, ::-1]
        Xrot, Yrot = np.meshgrid(y_edges, x_edges[::-1], indexing="ij")
        pcm = ax2d.pcolormesh(
            Xrot,
            Yrot,
            arr2d_rot,
            shading="auto",
            cmap="plasma",
            norm=norm,
            rasterized=True
        )
        ax2d.set_aspect("equal")

    if rays and ray_data is not None:
        plot_rays(
            ax2d,
            ray_data,
            xscale=ray_scale,
            yscale=ray_scale,
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

    cbar = fig.colorbar(
        pcm,
        ax=ax2d,
        orientation="horizontal",
        location="top",
        fraction=0.05,
        pad=0.12,
        aspect=30
    )
    cbar.set_label(field)

    # ============================================================
    # Bottom: 1D hydro + laser deposition lineout only
    # ============================================================
    c0, c1 = ["green", "red"]

    l1, = ax1d.plot(coord_plot, ni, lw=2, color=c0, label=r"$n_i$")
    l2, = ax1d.plot(coord_plot, ne, lw=2, color=c1, label=r"$n_e$")

    # Optional Te/Ti plotting block. Turn plot_temperature_lineout=True above.
    # if plot_temperature_lineout:
    #     tele_keV = tele * 1e-3
    #     tion_keV = tion * 1e-3
    #     l1.set_ydata(tele_keV)
    #     l2.set_ydata(tion_keV)
    #     l1.set_label(r"$T_e$ (keV)")
    #     l2.set_label(r"$T_i$ (keV)")
    #     ax1d.set_ylim(0, max(np.nanmax(tele_keV), np.nanmax(tion_keV)) * 1.2)

    finite_ni_ne = np.concatenate([ni[np.isfinite(ni)], ne[np.isfinite(ne)]])
    if finite_ni_ne.size > 0:
        ymin = 0.8 * np.nanmin(finite_ni_ne)
        ymax = 1.2 * np.nanmax(finite_ni_ne)

        # Avoid weird flat/singular y-limits.
        if ymin == ymax:
            ymin *= 0.9
            ymax *= 1.1

        ax1d.set_ylim(ymin, ymax)

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

    # Match 1D lineout x-limits to the horizontal axis of the rotated 2D image
    ax1d.set_xlim(ax2d.get_xlim())

    plt.tight_layout()

    # Get positions after tight_layout
    pos2d = ax2d.get_position()
    pos1d = ax1d.get_position()
    posc = cbar.ax.get_position()

    # Match the actual physical axis width of the 1D lineout to the 2D image
    ax1d.set_position([
        pos2d.x0,
        pos1d.y0,
        pos2d.width,
        pos1d.height
    ])

    # Since ax1d_2 is a twinx axis, move it too
    ax1d_2.set_position(ax1d.get_position())

    # Match colorbar width to the 2D image width too
    cbar.ax.set_position([
        pos2d.x0,
        posc.y0,
        pos2d.width,
        posc.height
    ])

    if savePlots:
        saveDir = Path(saveDir)
        saveDir.mkdir(parents=True, exist_ok=True)

        out = saveDir / f"{Path(fp).stem}_{field}_2D_with_lineout.png"
        plt.savefig(out, dpi=save_dpi, bbox_inches="tight")
        print(f"Saved {out}")

    plt.show()

    return fig, (ax2d, ax1d, ax1d_2)
