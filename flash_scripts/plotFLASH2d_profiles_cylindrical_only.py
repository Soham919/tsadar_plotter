import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from pathlib import Path
import yt
from flash_helpers import *
from scipy import constants

# Constants
mp = constants.proton_mass

# Atomic masses used for approximate number densities.
# Update these when you switch to real Si3N4 and H/He mixed EOS tables.
A_GAS = 1      # H gas placeholder
A_TARG = 28    # Si target placeholder
A_CHAM = 4     # He chamber placeholder


def _auto_downsample_2d(arr, max_pixels=1200):
    """
    Downsample a 2D array for display only.
    The lineout still uses the full-resolution arrays.
    """
    n0, n1 = arr.shape
    stride = max(1, int(np.ceil(max(n0, n1) / max_pixels)))
    return arr[::stride, ::stride], stride


def detect_first_density_jump_right_to_left(
    coord,
    dens,
    jump_threshold,
    window=5,
    mode="up",
    min_coord=None,
    max_coord=None,
):
    """
    Detect the first density jump while scanning from high coord to low coord.

    For cylindrical FLASH, this is usually used on a z-lineout near r=0,
    so it scans from large z to small z.
    """
    coord = np.asarray(coord)
    dens = np.asarray(dens)

    if coord.ndim != 1 or dens.ndim != 1:
        raise ValueError("coord and dens must be 1D arrays.")
    if len(coord) != len(dens):
        raise ValueError("coord and dens must have the same length.")

    good = np.isfinite(coord) & np.isfinite(dens)
    if min_coord is not None:
        good &= coord >= min_coord
    if max_coord is not None:
        good &= coord <= max_coord

    coord = coord[good]
    dens = dens[good]

    if len(coord) < 2 * window + 2:
        return {
            "found": False,
            "shock_coord": np.nan,
            "shock_index": None,
            "jump_value": np.nan,
            "dens_left": np.nan,
            "dens_right": np.nan,
        }

    # Sort so coord increases left-to-right numerically.
    sort_idx = np.argsort(coord)
    coord = coord[sort_idx]
    dens = dens[sort_idx]
    n = len(coord)

    # Scan from high coord to low coord.
    for i in range(n - window - 1, window - 1, -1):
        # Smaller-coordinate side
        dens_left = np.mean(dens[i - window + 1 : i + 1])
        # Larger-coordinate side
        dens_right = np.mean(dens[i + 1 : i + 1 + window])

        # Jump encountered when moving right-to-left.
        jump = dens_left - dens_right

        if mode == "up":
            condition = jump >= jump_threshold
        elif mode == "down":
            condition = jump <= -jump_threshold
        elif mode == "abs":
            condition = abs(jump) >= jump_threshold
        else:
            raise ValueError("mode must be 'up', 'down', or 'abs'.")

        if condition:
            shock_coord = 0.5 * (coord[i] + coord[i + 1])
            return {
                "found": True,
                "shock_coord": shock_coord,
                "shock_index": i,
                "jump_value": jump,
                "dens_left": dens_left,
                "dens_right": dens_right,
            }

    return {
        "found": False,
        "shock_coord": np.nan,
        "shock_index": None,
        "jump_value": np.nan,
        "dens_left": np.nan,
        "dens_right": np.nan,
    }


def plotFLASH2d_profiles_cylindrical(
    ds,
    ftype,
    field,
    useMicrons,
    savePlots,
    saveDir,
    fp,
    rays=False,
):
    """
    Cylindrical-only FLASH plotting function.

    Assumptions
    -----------
    FLASH/yt data axes are interpreted as:
        axis 0 -> r
        axis 1 -> z
        axis 2 -> theta, ignored here

    Top panel:
        2D r-z field image with z horizontal and r vertical.

    Bottom panel:
        1D lineout of ni, ne, and laser deposition.

    Useful default:
        lineout_axis = "z" and lineout_pos_um = 0 gives an axial lineout
        near the cylindrical axis, r ~= 0.
    """

    # ============================================================
    # Hardcoded settings
    # ============================================================
    lineout_axis = "z"       # "z" -> profile vs z at fixed r; "r" -> profile vs r at fixed z
    lineout_pos_um = 0       # fixed r or fixed z position, depending on lineout_axis

    max_plot_pixels = 1200
    plot_temperature_lineout = False
    plot_velocity_lineout = False
    save_dpi = 150

    # Shock detector settings
    detect_shock_position = True
    shock_density_source = "ne"   # "ne" or "ni"
    shock_jump_threshold = 5e19
    shock_window = 5
    shock_jump_mode = "up"        # "up", "down", or "abs"
    shock_min_coord = None
    shock_max_coord = None

    # Laser deposition quick check along the selected lineout coordinate.
    check_coord_um = 1000.0
    tol = 1e-30

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
    elif plot_velocity_lineout:
        velr2d = get_field_array(cg, "flash", "velx")[:, :, 0]  # FLASH velx -> radial velocity in cylindrical
        velz2d = get_field_array(cg, "flash", "vely")[:, :, 0]  # FLASH vely -> axial velocity in cylindrical

    print(f"{field} 2D shape:", arr2d.shape)
    print(f"{field} min/max:", np.nanmin(arr2d), np.nanmax(arr2d))

    # ============================================================
    # Coordinates: FLASH axes are interpreted as r, z, theta
    # ============================================================
    r0, r1 = float(ds.domain_left_edge[0]), float(ds.domain_right_edge[0])
    z0, z1 = float(ds.domain_left_edge[1]), float(ds.domain_right_edge[1])

    nr, nz = arr2d.shape

    r_edges_cm = np.linspace(r0, r1, nr + 1)
    z_edges_cm = np.linspace(z0, z1, nz + 1)

    r_centers_cm = 0.5 * (r_edges_cm[:-1] + r_edges_cm[1:])
    z_centers_cm = 0.5 * (z_edges_cm[:-1] + z_edges_cm[1:])

    if useMicrons:
        scale = 1e4
        r_edges = r_edges_cm * scale
        z_edges = z_edges_cm * scale
        r_centers = r_centers_cm * scale
        z_centers = z_centers_cm * scale
        lineout_pos = lineout_pos_um
        xlabel_r = r"$r$ [$\mu$m]"
        xlabel_z = r"$z$ [$\mu$m]"
        coord_unit = r"$\mu$m"
        ray_scale = scale
    else:
        scale = 1.0
        r_edges = r_edges_cm
        z_edges = z_edges_cm
        r_centers = r_centers_cm
        z_centers = z_centers_cm
        lineout_pos = lineout_pos_um * 1e-4
        xlabel_r = r"$r$ [cm]"
        xlabel_z = r"$z$ [cm]"
        coord_unit = "cm"
        ray_scale = scale

    # ============================================================
    # Extract 1D lineout
    # ============================================================
    laxis = lineout_axis.lower()

    if laxis == "z":
        # Axial lineout at fixed radius. For lineout_pos_um=0, this picks the
        # first cell center near the cylindrical axis, not exactly r=0.
        ir = np.argmin(np.abs(r_centers - lineout_pos))

        coord_plot = z_centers
        xlabel_line = xlabel_z

        dens = dens2d[ir, :]
        depo = depo2d[ir, :]
        ye   = ye2d[ir, :]

        gas  = gas2d[ir, :]
        targ = targ2d[ir, :]
        cham = cham2d[ir, :]

        if plot_temperature_lineout:
            tele = tele2d[ir, :]
            tion = tion2d[ir, :]
        elif plot_velocity_lineout:
            velr = velr2d[ir, :]
            velz = velz2d[ir, :]

        line_location = r_centers[ir]
        line_title = rf"Lineout vs $z$ at $r = {line_location:.3f}$ {coord_unit}"
        fixed_coord_name = "r"
        lineout_coord_name = "z"

    elif laxis == "r":
        # Radial lineout at fixed axial position.
        iz = np.argmin(np.abs(z_centers - lineout_pos))

        coord_plot = r_centers
        xlabel_line = xlabel_r

        dens = dens2d[:, iz]
        depo = depo2d[:, iz]
        ye   = ye2d[:, iz]

        gas  = gas2d[:, iz]
        targ = targ2d[:, iz]
        cham = cham2d[:, iz]

        if plot_temperature_lineout:
            tele = tele2d[:, iz]
            tion = tion2d[:, iz]
        elif plot_velocity_lineout:
            velr = velr2d[:, iz]
            velz = velz2d[:, iz]

        line_location = z_centers[iz]
        line_title = rf"Lineout vs $r$ at $z = {line_location:.3f}$ {coord_unit}"
        fixed_coord_name = "z"
        lineout_coord_name = "r"

    else:
        raise ValueError("For this cylindrical script, lineout_axis must be 'r' or 'z'.")

    print(f"Cylindrical lineout axis: {lineout_axis}")
    print(f"Requested fixed {fixed_coord_name} position: {lineout_pos_um} microns")
    print(f"Nearest fixed {fixed_coord_name} location used: {line_location:.4f} {coord_unit}")

    # ============================================================
    # Convert lineout fields
    # ============================================================
    mp_g = mp * 1000.0

    gas_dens    = dens * gas  / (A_GAS * mp_g)
    target_dens = dens * targ / (A_TARG * mp_g)
    cham_dens   = dens * cham / (A_CHAM * mp_g)

    ni = gas_dens + target_dens + cham_dens
    ne = ye * dens / mp_g
    las_depo = dens * depo * 1e-7

    # ============================================================
    # Detect shock position from 1D lineout
    # ============================================================
    shock_result = None
    shock_coord = np.nan

    if detect_shock_position:
        if shock_density_source.lower() == "ne":
            shock_dens = ne
        elif shock_density_source.lower() == "ni":
            shock_dens = ni
        else:
            raise ValueError("shock_density_source must be either 'ne' or 'ni'.")

        shock_result = detect_first_density_jump_right_to_left(
            coord_plot,
            shock_dens,
            jump_threshold=shock_jump_threshold,
            window=shock_window,
            mode=shock_jump_mode,
            min_coord=shock_min_coord,
            max_coord=shock_max_coord,
        )
        shock_coord = shock_result["shock_coord"]

        if shock_result["found"]:
            print(
                f"Detected shock using {shock_density_source}: "
                f"{shock_coord:.4f} {coord_unit}, "
                f"jump = {shock_result['jump_value']:.3e}, "
                f"t = {sim_time_ns:.4f} ns"
            )
        else:
            print(
                f"No shock detected using {shock_density_source} "
                f"at t = {sim_time_ns:.4f} ns. "
                f"Try lowering shock_jump_threshold."
            )

    # ============================================================
    # Quick check: is laser deposition nonzero near chosen lineout coordinate?
    # ============================================================
    RED = "\033[91m"
    RESET = "\033[0m"

    i_check = np.argmin(np.abs(coord_plot - check_coord_um))
    val_check = las_depo[i_check]

    if abs(val_check) > tol:
        print(
            RED
            + f"Laser deposition is NONZERO near {lineout_coord_name} = {check_coord_um:.3g} um "
            + f"(nearest {lineout_coord_name} = {coord_plot[i_check]:.3g} um, "
            + f"value = {val_check:.3e} J/cm^3)"
            + RESET
        )
    else:
        print(
            RED
            + f"Laser deposition is zero near {lineout_coord_name} = {check_coord_um:.3g} um "
            + f"(nearest {lineout_coord_name} = {coord_plot[i_check]:.3g} um, "
            + f"value = {val_check:.3e} J/cm^3)"
            + RESET
        )

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
        2,
        1,
        figsize=(9, 7),
        gridspec_kw={"height_ratios": [2.2, 1]},
        constrained_layout=False,
    )

    # ============================================================
    # Top: 2D r-z plot
    # ============================================================
    positive = arr2d[arr2d > 0]
    if positive.size > 0:
        norm = LogNorm(vmin=np.nanmin(positive), vmax=np.nanmax(arr2d))
    else:
        print(f"Warning: {field} has no positive values. Using linear scale.")
        norm = Normalize(vmin=np.nanmin(arr2d), vmax=np.nanmax(arr2d))

    # Display z horizontally and r vertically.
    # arr2d has shape [r, z]. imshow expects [vertical, horizontal].
    # Use arr2d[::-1, :] so r increases upward on the displayed axis.
    arr_img = arr2d[::-1, :]
    arr_img_plot, display_stride = _auto_downsample_2d(arr_img, max_pixels=max_plot_pixels)

    if display_stride > 1:
        print(
            f"Display downsampling enabled: showing every {display_stride}th cell "
            f"for 2D image only. Lineout remains full resolution."
        )

    pcm = ax2d.imshow(
        arr_img_plot,
        extent=[z_edges[0], z_edges[-1], r_edges[0], r_edges[-1]],
        origin="upper",
        aspect="equal",
        interpolation="nearest",
        cmap="inferno",
        norm=norm,
    )

    if rays and ray_data is not None:
        # Assumes ray data are stored as x/y-like coordinates corresponding to r/z.
        # If they look swapped, swap xscale/yscale handling inside plot_rays.
        plot_rays(ax2d, ray_data, xscale=ray_scale, yscale=ray_scale, color="w", lw=0.8)

    # Draw selected lineout location on 2D r-z image.
    if laxis == "z":
        # Fixed r -> horizontal line on r-z image.
        ax2d.axhline(line_location, color="w", linestyle="--", linewidth=1.2)
    else:
        # Fixed z -> vertical line on r-z image.
        ax2d.axvline(line_location, color="w", linestyle="--", linewidth=1.2)

    ax2d.set_xlabel(xlabel_z)
    ax2d.set_ylabel(xlabel_r)
    ax2d.set_title(f"{field}, t = {sim_time_ns:.3f} ns")

    cbar = fig.colorbar(
        pcm,
        ax=ax2d,
        orientation="horizontal",
        location="top",
        fraction=0.05,
        pad=0.12,
        aspect=30,
    )
    cbar.set_label(field)

    # ============================================================
    # Bottom: 1D lineout
    # ============================================================
    l1, = ax1d.plot(coord_plot, ni, lw=2, color="green", label=r"$n_i$")
    l2, = ax1d.plot(coord_plot, ne, lw=2, color="red", label=r"$n_e$")

    if detect_shock_position and shock_result is not None and shock_result["found"]:
        ax1d.axvline(
            shock_coord,
            color="k",
            linestyle=":",
            linewidth=2,
            label=rf"shock = {shock_coord:.1f} {coord_unit}",
        )

    if plot_temperature_lineout:
        tele_keV = tele * 1e-3
        tion_keV = tion * 1e-3
        l1.set_ydata(tele_keV)
        l2.set_ydata(tion_keV)
        l1.set_label(r"$T_e$ (keV)")
        l2.set_label(r"$T_i$ (keV)")
        ax1d.set_ylim(0, max(np.nanmax(tele_keV), np.nanmax(tion_keV)) * 1.2)
        ax1d.set_ylabel(r"$T$ [keV]")

    elif plot_velocity_lineout:
        l1.set_ydata(velr)
        l2.set_ydata(velz)
        l1.set_label(r"$v_r$ (cm/s)")
        l2.set_label(r"$v_z$ (cm/s)")
        ax1d.set_ylim(0, max(np.nanmax(velr), np.nanmax(velz)) * 1.2)
        ax1d.set_ylabel(r"$v$ [cm/s]")

    else:
        finite_ni_ne = np.concatenate([ni[np.isfinite(ni)], ne[np.isfinite(ne)]])
        if finite_ni_ne.size > 0:
            ymin = 0.8 * np.nanmin(finite_ni_ne)
            ymax = 1.2 * np.nanmax(finite_ni_ne)
            if ymin == ymax:
                ymin *= 0.9
                ymax *= 1.1
            ax1d.set_ylim(ymin, ymax)
        ax1d.set_ylabel(r"$n$ [cm$^{-3}$]")

    ax1d.set_xlabel(xlabel_line)
    ax1d.set_title(line_title)
    ax1d.grid(True, alpha=0.3)

    ax1d_2 = ax1d.twinx()
    l3, = ax1d_2.plot(
        coord_plot[50:-1],
        las_depo[50:-1],
        color="blueviolet",
        linewidth=2,
        linestyle="--",
        label=r"$E_{dep}$",
    )

    ax1d_2.set_ylabel(r"$E_{dep}$ [J/cm$^3$]", color="blueviolet")
    ax1d_2.tick_params(axis="y", colors="blueviolet")
    ax1d_2.spines["right"].set_color("blueviolet")
    ax1d_2.spines["right"].set_linewidth(1.5)

    lines = [l1, l2, l3]
    labels = [line.get_label() for line in lines]
    ax1d.legend(lines, labels, loc="upper right", frameon=False)

    # Match 1D lineout x-limits to the corresponding horizontal 2D axis only for z-lineouts.
    if laxis == "z":
        ax1d.set_xlim(ax2d.get_xlim())
    else:
        ax1d.set_xlim(coord_plot[0], coord_plot[-1])

    plt.tight_layout()

    # Align 1D axis and colorbar width with the 2D axis width.
    pos2d = ax2d.get_position()
    pos1d = ax1d.get_position()
    posc = cbar.ax.get_position()

    ax1d.set_position([pos2d.x0, pos1d.y0, pos2d.width, pos1d.height])
    ax1d_2.set_position(ax1d.get_position())
    cbar.ax.set_position([pos2d.x0, posc.y0, pos2d.width, posc.height])

    if savePlots:
        saveDir = Path(saveDir)
        saveDir.mkdir(parents=True, exist_ok=True)
        out = saveDir / f"{Path(fp).stem}_{field}_cyl_rz_with_lineout.png"
        plt.savefig(out, dpi=save_dpi, bbox_inches="tight")
        print(f"Saved {out}")

    return fig, (ax2d, ax1d, ax1d_2), shock_coord, shock_result
