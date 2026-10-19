import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from pathlib import Path
import sys
import yt
import glob
import os
import re
from scipy import constants
from scipy.ndimage import gaussian_filter
# ----- Mac -----#
#sys.path.append("/Users/soham/Documents/Kinshock/Scripts")
# ----- Windows -----#
sys.path.append(r"\\profiles\Users$\sban\Documents\Scripts")
from plasma_tools import mfp_ii
from shock_detector_1d import detect_simple_density_jump

#--------- Constants ---------- #
mp = constants.proton_mass # mass of H(kg)
eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi
c = constants.c

# ============================================================
# Helper functions
# ============================================================

def get_flash_sim_time(filename):
    """
    Read FLASH simulation time directly from HDF5 metadata.
    Returns time in seconds.
    """
    with h5py.File(filename, "r") as f:
        real_scalars = f["real scalars"][:]

        for row in real_scalars:
            name = row[0]
            value = row[1]

            if isinstance(name, bytes):
                name = name.decode("utf-8").strip()
            else:
                name = str(name).strip()

            if name == "time":
                return float(value)

    raise ValueError(f"No time entry found in {filename}")


def get_file_number(filename):
    """
    Extract number from something like:
    ks_hdf5_plt_cnt_0042
    """
    match = re.search(r"plt_cnt_(\d+)", os.path.basename(filename))
    if match:
        return int(match.group(1))
    return -1


def make_flash_time_table(
    fp,
    file_pattern="*plt_cnt_*",
    output_txt="flash_plotfile_times.txt"
):
    """
    Scan FLASH plot files once and save file number/time mapping.
    Does NOT use os.chdir().
    """

    import os
    import glob

    search_pattern = os.path.join(fp, file_pattern)
    output_path = os.path.join(fp, output_txt)

    files = sorted(
        glob.glob(search_pattern),
        key=get_file_number
    )

    if len(files) == 0:
        raise FileNotFoundError(f"No files found matching pattern: {search_pattern}")

    rows = []

    for file in files:
        try:
            file_num = get_file_number(file)
            sim_time_s = get_flash_sim_time(file)
            sim_time_ns = sim_time_s * 1e9

            # Save only the filename, not the full path
            filename_only = os.path.basename(file)

            rows.append((file_num, sim_time_s, sim_time_ns, filename_only))

            print(
                f"{file_num:05d}  "
                f"{sim_time_s:.12e} s  "
                f"{sim_time_ns:.6f} ns  "
                f"{filename_only}"
            )

        except Exception as e:
            print(f"Skipping {file}: {e}")

    if len(rows) == 0:
        raise RuntimeError("No valid FLASH plot files found.")

    with open(output_path, "w") as f:
        f.write("# file_number sim_time_s sim_time_ns filename\n")

        for file_num, sim_time_s, sim_time_ns, filename in rows:
            f.write(
                f"{file_num:d} "
                f"{sim_time_s:.12e} "
                f"{sim_time_ns:.12e} "
                f"{filename}\n"
            )

    print(f"\nSaved time table to: {output_path}")

    return rows


def find_flash_file_from_number(file_number, fp):
    """
    Find FLASH plotfile(s) based on the final 4-digit file number.

    Parameters
    ----------
    file_number : int or array-like
        File number(s), e.g. 2 or [2, 10, 25].

    fp : str or Path
        Directory containing FLASH plotfiles.

    Returns
    -------
    Path
        If a single file number is supplied.

    list of Path
        If multiple file numbers are supplied.
    """

    fp = Path(fp)

    # Was the input a single number?
    input_is_scalar = np.isscalar(file_number)

    file_numbers = np.atleast_1d(file_number)

    results = []

    for num in file_numbers:

        num = int(num)

        # FLASH convention: _hdf5_plt_cnt_XXXX
        pattern = f"*hdf5_plt_cnt_{num:04d}"

        matches = list(fp.glob(pattern))

        if len(matches) == 0:
            raise FileNotFoundError(
                f"Could not find FLASH plotfile number {num} "
                f"in {fp}"
            )

       
        results.extend(matches)

    # Single number -> single Path
    if input_is_scalar:
        return results[0]

    # Array -> list of Paths
    return results


def find_nearest_flash_file_from_table(
    target_time_ns,
    fp,
    time_table="flash_plotfile_times.txt",
):
    """
    Find nearest FLASH plot file(s) using precomputed time table.

    target_time_ns can be:
      - a single number, like 1.0
      - a list/array, like [0.5, 1.0, 1.5]
    """

    import os
    import numpy as np

    table_path = os.path.join(fp, time_table)

    file_numbers = []
    sim_times_s = []
    sim_times_ns = []
    filenames = []

    with open(table_path, "r") as f:
        for line in f:
            if line.strip().startswith("#") or len(line.strip()) == 0:
                continue

            parts = line.split()

            file_numbers.append(int(parts[0]))
            sim_times_s.append(float(parts[1]))
            sim_times_ns.append(float(parts[2]))

            # Return full path to file
            filenames.append(os.path.join(fp, parts[3]))

    sim_times_s = np.array(sim_times_s)
    sim_times_ns = np.array(sim_times_ns)
    file_numbers = np.array(file_numbers)

    # Detect whether input is scalar or array-like
    input_is_scalar = np.isscalar(target_time_ns)
    target_times_ns = np.atleast_1d(target_time_ns)

    results = []

    for t_ns in target_times_ns:
        idx = np.argmin(np.abs(sim_times_ns - t_ns))

        nearest_file = filenames[idx]
        nearest_file_number = file_numbers[idx]
        nearest_time_s = sim_times_s[idx]
        nearest_time_ns = sim_times_ns[idx]

        print("Requested time:")
        print(f"  {t_ns:.6f} ns")

        print("Nearest FLASH plot file:")
        print(f"  file name   = {nearest_file}")
        print(f"  file number = {nearest_file_number}")
        print(f"  sim time    = {nearest_time_ns:.6f} ns")
        print(f"  difference  = {abs(nearest_time_ns - t_ns):.6f} ns")
        print()

        results.append(
            (nearest_file, nearest_file_number, nearest_time_s, nearest_time_ns)
        )

    # Preserve old behavior for a single input
    if input_is_scalar:
        return results[0]

    return results


def get_sim_time_ns(ds):
    """
    FLASH time is usually in seconds for these setups.
    """
    sim_time = float(ds.current_time)
    return sim_time * 1e9


def print_fields(ds, ad):
    print("\nNative fields:")
    for field in ds.field_list:
        try:
            print("  ", f"{field}: {ad[field].units}")
        except Exception as e:
            print(f"  {field}: failed ({e})")

    print("\nDerived fields:")
    for field in ds.derived_field_list:
        try:
            print("  ", f"{field}: {ad[field].units}")
        except Exception as e:
            print(f"  {field}: failed ({type(e).__name__}: {e})")


def get_covering_grid(ds):
    """
    Returns max-level covering grid and dimensions,
    avoiding fake FLASH dimensions for 1D/2D runs.
    """
    level = ds.index.max_level
    base_dims = np.array(ds.domain_dimensions, dtype=int)
    refined_dims = base_dims * (2 ** level)

    plot_dim = detect_plot_dim_from_ds(ds)

    dims = np.ones(3, dtype=int)

    if plot_dim == "1d":
        dims[0] = refined_dims[0]

    elif plot_dim == "2d":
        dims[0] = refined_dims[0]
        dims[1] = refined_dims[1]

    elif plot_dim == "3d":
        dims = refined_dims

    else:
        raise ValueError(f"Unknown plot_dim: {plot_dim}")

    print("max level:", level)
    print("domain_dimensions:", ds.domain_dimensions)
    print("detected plot_dim:", plot_dim)
    print("raw refined dims:", refined_dims)
    print("safe covering grid dims:", dims)

    cg = ds.covering_grid(
        level=level,
        left_edge=ds.domain_left_edge,
        dims=dims
    )

    return cg, dims


def detect_plot_dim_from_ds(ds):

    """

    Detect actual FLASH dimensionality from base domain_dimensions,

    not covering-grid dims.

    yt covering_grid multiplies even singleton dimensions by 2**level,

    so a true 1D run like [2048, 1, 1] can become [16384, 8, 8].

    """

    base_dims = np.array(ds.domain_dimensions)

    active_dims = np.sum(base_dims > 1)

    print("base active dimensions:", active_dims)

    if active_dims <= 1:

        return "1d"

    elif active_dims == 2:

        return "2d"

    else:

        return "3d"


def coord_to_plot_units(x_cm, useMicrons):
    """
    FLASH uses cgs lengths, so coordinates are usually cm.
    Convert cm to microns if requested.
    """
    if useMicrons:
        return x_cm * 1e4, r"$x$ [$\mu$m]"
    else:
        return x_cm, r"$x$ [cm]"


def get_field_array(cg, ftype, field):
    """
    Safely get a FLASH field from covering grid.
    """
    return cg[(ftype, field)].to_ndarray()


def squeeze_to_1d(arr):

    """

    Convert FLASH/yt covering-grid array to a 1D profile.

    For true 1D FLASH, yt may still return shape:

        (nx, 8, 8)

    after max-level covering_grid because singleton dims got refined.

    We take the first y,z pencil:

        arr[:, 0, 0]

    """

    arr = np.asarray(arr)

    if arr.ndim == 1:

        return arr

    if arr.ndim == 3:

        return arr[:, 0, 0]

    if arr.ndim == 2:

        return arr[:, 0]

    raise ValueError(f"Cannot convert array with shape {arr.shape} to 1D.")


def detect_first_density_jump_right_to_left(
    coord, sim_time, dens, tele, tion, vely,
    gas, targ, cham,
    jump_threshold,
    window=5,
    post_width=1000,
    post_offset=0,
    gas_threshold=0.90,
    min_post_cells=1,
    mode="up",
    min_coord=None,
    max_coord=None,
):
    """
    Detect the first gas-density jump scanning right to left.

    Then find the maximum density in the contiguous gas-rich
    region immediately behind the shock.

    All post-shock properties are extracted at that same cell.

    post_width and post_offset use the same units as coord.

    statistic and sampling_mode are retained for compatibility
    with existing calls but are no longer used.
    """

    # --------------------------------------------------------
    # Input arrays and validation
    # --------------------------------------------------------
    arrays = [
        np.asarray(a, dtype=float)
        for a in (coord, dens, tele, tion, vely, gas, targ, cham)
    ]

    coord, dens, tele, tion, vely, gas, targ, cham = arrays

    if any(a.ndim != 1 for a in arrays):
        raise ValueError("All input arrays must be 1D.")

    if any(len(a) != len(coord) for a in arrays):
        raise ValueError("All input arrays must have equal length.")

    if window < 1:
        raise ValueError("window must be >= 1.")

    if post_width <= post_offset or post_offset < 0:
        raise ValueError("Require 0 <= post_offset < post_width.")

    if min_post_cells < 1:
        raise ValueError("min_post_cells must be >= 1.")

    if mode not in ("up", "down", "abs"):
        raise ValueError("mode must be 'up', 'down', or 'abs'.")

    # --------------------------------------------------------
    # Filter invalid values
    # --------------------------------------------------------
    good = np.logical_and.reduce([
        np.isfinite(a) for a in arrays
    ])

    if min_coord is not None:
        good &= coord >= min_coord

    if max_coord is not None:
        good &= coord <= max_coord

    coord, dens, tele, tion, vely, gas, targ, cham = (
        a[good] for a in arrays
    )

    # --------------------------------------------------------
    # Default result: preserve all original return variables
    # --------------------------------------------------------
    result = {
        "found": False,
        "shock_coord": np.nan,
        "Time": sim_time,
        "shock_index": None,
        "jump_value": np.nan,
        "dens_left": np.nan,
        "dens_right": np.nan,
        "n_post_cells": 0,
        "postshock_dens": np.nan,
        "postshock_tele": np.nan,
        "postshock_tion": np.nan,
        "postshock_vely": np.nan,
        "postshock_gas_fraction": np.nan,
        "postshock_target_fraction": np.nan,
        "postshock_chamber_fraction": np.nan,
        "postshock_coord_min": np.nan,
        "postshock_coord_max": np.nan,
        "postshock_peak_coord": np.nan,
    }

    if len(coord) < 2 * window + 2:
        return result

    # --------------------------------------------------------
    # Sort coordinates low -> high
    # --------------------------------------------------------
    order = np.argsort(coord)

    coord, dens, tele, tion, vely, gas, targ, cham = (
        a[order]
        for a in (coord, dens, tele, tion, vely, gas, targ, cham)
    )

    n = len(coord)

    # --------------------------------------------------------
    # Scan right to left: original shock detection
    # --------------------------------------------------------
    for i in range(n - window - 1, window - 1, -1):

        left_slice = slice(i - window + 1, i + 1)
        right_slice = slice(i + 1, i + 1 + window)

        # Both detection windows must be gas-rich.
        if (
            np.any(gas[left_slice] < gas_threshold)
            or np.any(gas[right_slice] < gas_threshold)
        ):
            continue

        dens_left = np.mean(dens[left_slice])
        dens_right = np.mean(dens[right_slice])

        jump = dens_left - dens_right

        if mode == "up":
            condition = jump >= jump_threshold
        elif mode == "down":
            condition = jump <= -jump_threshold
        else:
            condition = abs(jump) >= jump_threshold

        if not condition:
            continue

        # ----------------------------------------------------
        # Shock detected
        # ----------------------------------------------------
        shock_coord = 0.5 * (coord[i] + coord[i + 1])

        result.update({
            "found": True,
            "shock_coord": float(shock_coord),
            "shock_index": int(i),
            "jump_value": float(jump),
            "dens_left": float(dens_left),
            "dens_right": float(dens_right),
        })

        # ----------------------------------------------------
        # Find contiguous gas-rich region behind shock
        # ----------------------------------------------------
        post_indices = []

        for j in range(i, -1, -1):

            distance = shock_coord - coord[j]

            if distance > post_width:
                break

            # Stop when we reach target/mixed material.
            if gas[j] < gas_threshold:
                break

            # Skip shock ramp if requested.
            if distance < post_offset:
                continue

            post_indices.append(j)

        result["n_post_cells"] = len(post_indices)

        if len(post_indices) < min_post_cells:
            return result

        post_indices = np.asarray(post_indices, dtype=int)

        # ----------------------------------------------------
        # Find maximum density behind shock
        # ----------------------------------------------------
        peak_idx = post_indices[
            np.argmax(dens[post_indices])
        ]

        # ----------------------------------------------------
        # Extract all properties at the density peak
        # ----------------------------------------------------
        result.update({
            "postshock_dens": float(dens[peak_idx]),
            "postshock_tele": float(tele[peak_idx]),
            "postshock_tion": float(tion[peak_idx]),
            "postshock_vely": float(vely[peak_idx]),
            "postshock_gas_fraction": float(gas[peak_idx]),
            "postshock_target_fraction": float(targ[peak_idx]),
            "postshock_chamber_fraction": float(cham[peak_idx]),
            "postshock_coord_min": float(
                np.min(coord[post_indices])
            ),
            "postshock_coord_max": float(
                np.max(coord[post_indices])
            ),
            "postshock_peak_coord": float(coord[peak_idx]),
        })

        return result

    return result

def plot_1d_field(
    ds,
    field,
    file_number,
    savePlots=False,
    saveDir=None,
    shock_jump_threshold=5.0e-05,
    post_offset=0.0,
    detect_shock=True,
):
    """
    Plot a single 1D FLASH field and optionally detect
    the shock using the density profile.

    post_offset is in cm.
    """

    # --------------------------------------------------
    # Load FLASH data
    # --------------------------------------------------
    t = get_sim_time_ns(ds)
    cg, dims = get_covering_grid(ds)

    shock_result = None

    # Check requested field
    if (
        ("flash", field) not in ds.field_list
        and ("flash", field) not in ds.derived_field_list
    ):
        raise ValueError(f"Field '{field}' not found in dataset.")

    # Field to plot
    y = squeeze_to_1d(
        get_field_array(cg, "flash", field)
    )

    # --------------------------------------------------
    # Radial coordinate
    # --------------------------------------------------
    nx = y.size

    x0 = float(ds.domain_left_edge[0])
    x1 = float(ds.domain_right_edge[0])

    dx = (x1 - x0) / nx

    x = x0 + (np.arange(nx) + 0.5) * dx

    # --------------------------------------------------
    # Shock detection (independent of plotted field)
    # --------------------------------------------------
    if detect_shock:

        dens = squeeze_to_1d(
            get_field_array(cg, "flash", "dens")
        )

        tele = squeeze_to_1d(
            get_field_array(cg, "flash", "tele")
        )

        tion = squeeze_to_1d(
            get_field_array(cg, "flash", "tion")
        )

        velx = squeeze_to_1d(
            get_field_array(cg, "flash", "velx")
        )

        shock_result = detect_simple_density_jump(
            coord=x,
            sim_time=t,
            dens=dens,
            tele=tele,
            tion=tion,
            velx=velx,
            jump_threshold=shock_jump_threshold,
            post_offset=post_offset,
        )

    # --------------------------------------------------
    # Plot requested field
    # --------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 4.5))

    ax.plot(x, y, linewidth=2)

    if shock_result is not None and shock_result["found"]:

        ax.axvline(
            shock_result["shock_coord"]*10**(-4),
            color="red",
            linestyle="--",
            linewidth=2,
            label="Detected Shock"
        )

        ax.legend()
    ax.set_yscale("log")
    ax.set_xlabel("x [cm]")
    ax.set_ylabel(field)
    ax.set_title(f"{field}, file = {file_number}")

    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    #plt.show()

    # --------------------------------------------------
    # Save figure
    # --------------------------------------------------
    if savePlots:

        saveDir = Path(saveDir)
        saveDir.mkdir(parents=True, exist_ok=True)

        number = str(file_number).split("_")[-1]

        output_path = saveDir / f"{number}_{field}_1D.png"

        fig.savefig(
            output_path,
            dpi=150,
            bbox_inches="tight"
        )

        print(f"Saved {output_path}")

    return fig, ax, t, shock_result



def plot_1d_profiles_with_material_fractions(ds, fields, useMicrons, savePlots, saveDir, fp, rays=False,):
    """
    Plot 1D FLASH profiles with material-fraction demarcation.
    """

    cg, dims = get_covering_grid(ds)
    sim_time_ns = get_sim_time_ns(ds)

    x0 = float(ds.domain_left_edge[0])
    x1 = float(ds.domain_right_edge[0])

    # ----------------------------------------------------------
    # Find a valid field to determine nx
    # ----------------------------------------------------------
    example_field = None

    for f in fields:
        if (("flash", f) in ds.field_list or ("flash", f) in ds.derived_field_list):
            example_field = f
            break

    if example_field is None:
        raise ValueError("None of the requested 1D fields were found in the dataset.")

    example_arr = squeeze_to_1d(get_field_array(cg, 'flash',example_field))
    nx = example_arr.size

    # ----------------------------------------------------------
    # Cell-center coordinates
    # ----------------------------------------------------------
    dx = (x1 - x0) / nx
    x_cm = x0 + (np.arange(nx) + 0.5) * dx
    x_plot, xlabel = coord_to_plot_units(x_cm, useMicrons,)

    # ----------------------------------------------------------
    # Material fractions
    # ----------------------------------------------------------
    chamber_line = squeeze_to_1d(get_field_array(cg, 'flash', "cham"))
    gas_line = squeeze_to_1d(get_field_array(cg, 'flash', "gas"))
    target_line = squeeze_to_1d(get_field_array(cg, 'flash', "targ"))

    # Make sure species arrays match hydro grid
    expected_shape = (nx,)

    for name, species_array in {"cham": chamber_line,"gas": gas_line,"targ": target_line,}.items():

        if species_array.shape != expected_shape:
            raise ValueError(f"{name} has shape {species_array.shape}, "f"but expected {expected_shape}.")

    # ----------------------------------------------------------
    # Loop over requested fields
    # ----------------------------------------------------------
    for field in fields:

        if (("flash", field) not in ds.field_list and ("flash", field) not in ds.derived_field_list):
            print(f"Skipping {field}: not found")
            continue

        y = squeeze_to_1d(get_field_array(cg, 'flash', field))

        if y.shape != expected_shape:
            raise ValueError(f"{field} has shape {y.shape}, "f"but expected {expected_shape}.")

        # ======================================================
        # Figure
        # ======================================================
        fig, (ax_field, ax_fraction) = plt.subplots(2, 1, figsize=(8, 6), sharex=True, gridspec_kw={"height_ratios": [2.5, 1.0],}, constrained_layout=True,)

        # ======================================================
        # Panel 1: FLASH field
        # ======================================================
        ax_field.plot(x_plot, y, linewidth=2,)
        ax_field.set_ylabel(field)
        ax_field.set_title(f"{field}, t = {sim_time_ns:.3f} ns")
        ax_field.grid(True, alpha=0.3,)

        ax_field.tick_params(axis="x",labelbottom=False,)

        ax_field.ticklabel_format(axis="y",style="sci",scilimits=(0, 0),)

        # ======================================================
        # Panel 2: material fractions
        # ======================================================
        ax_fraction.plot(x_plot,chamber_line,linewidth=2,label="Chamber",)

        ax_fraction.plot(x_plot,gas_line,linewidth=2,label="Gas",)

        ax_fraction.plot(x_plot,target_line,linewidth=2,label="Target",)

        ax_fraction.set_xlabel(xlabel)
        ax_fraction.set_ylabel("Fraction")
        ax_fraction.set_title("Material fractions")

        ax_fraction.set_ylim(-0.05,1.05,)

        ax_fraction.grid(True,alpha=0.3,)

        ax_fraction.legend(loc="best",frameon=False,)

        # ------------------------------------------------------
        # Match x ranges
        # ------------------------------------------------------
        xmin_plot = x_plot[0]
        xmax_plot = x_plot[-1]

        ax_field.set_xlim( xmin_plot, xmax_plot,)

        ax_fraction.set_xlim(xmin_plot,xmax_plot,)

        # ======================================================
        # Save
        # ======================================================
        if savePlots:

            saveDir = Path(saveDir)

            saveDir.mkdir(
                parents=True,
                exist_ok=True,
            )

            output_path = (
                saveDir
                / f"{Path(fp).stem}_{field}_1D.png"
            )

            fig.savefig(
                output_path,
                dpi=150,
                bbox_inches="tight",
            )

            print(f"Saved {output_path}")

       # plt.show()


def _auto_downsample_2d(arr, max_pixels=1200):
    """
    Downsample a 2D array for display only.
    The full-resolution array remains available for lineouts.
    """
    n0, n1 = arr.shape
    stride = max(1, int(np.ceil(max(n0, n1) / max_pixels)),)
    return arr[::stride, ::stride], stride


def plot_2d_profiles_with_material_fractions(
    ds,
    fp,
    ftype,
    field,
    lineout_dir="x",
    lineout_pos=0.0,
    useMicrons=True,
    savePlot=True,
    saveDir=Path("."),
    max_plot_pixels=1200,
):
    """
    Plot one FLASH field, a full-resolution lineout, and material fractions.

    Parameters
    ----------
    ds
        Loaded yt FLASH dataset.

    fp
        Path to the FLASH plotfile.

    ftype
        FLASH field type passed to get_field_array().

    field
        FLASH field name to plot.

    lineout_dir : {"x", "y"}
        Coordinate held fixed.

        lineout_dir="x":
            Hold x fixed at lineout_pos and plot versus y.

        lineout_dir="y":
            Hold y fixed at lineout_pos and plot versus x.

    lineout_pos : float
        Requested lineout position.

        In microns when useMicrons=True.
        In cm when useMicrons=False.

    useMicrons : bool
        Convert spatial coordinates from cm to microns.

    savePlot : bool
        Save the figure.

    saveDir
        Directory where the figure is saved.

    max_plot_pixels : int
        Maximum approximate pixel count along the longest dimension
        of the displayed 2D image.

        This affects only the image. Lineouts remain full resolution.

    Returns
    -------
    fig, axes, line_data
        axes = (ax_line, ax_2d, ax_fraction)

        line_data is a dictionary containing the extracted lineouts.
    """

    lineout_dir = lineout_dir.lower()

    if lineout_dir not in {"x", "y"}:
        raise ValueError(
            "lineout_dir must be either 'x' or 'y'."
        )

    # ==========================================================
    # Load covering grid and requested field
    # ==========================================================
    cg, dims = get_covering_grid(ds)
    sim_time_ns = get_sim_time_ns(ds)
    
    arr = get_field_array(cg, ftype, field,)
    arr = arr* kb / e if field in {"tele", "tion"} else arr
    arr2d = np.asarray(arr[:, :, 0])

    if arr2d.ndim != 2:
        raise ValueError(f"Expected a 2D field array, got shape {arr2d.shape}.")

    nx, ny = arr2d.shape

    print(f"{field} shape:", arr2d.shape)

    # ==========================================================
    # Coordinate arrays
    #
    # FLASH spatial coordinates are assumed to be in cm.
    # ==========================================================
    x0_cm = float(ds.domain_left_edge[0])
    x1_cm = float(ds.domain_right_edge[0])

    y0_cm = float(ds.domain_left_edge[1])
    y1_cm = float(ds.domain_right_edge[1])

    x_edges_cm = np.linspace( x0_cm, x1_cm, nx + 1,)

    y_edges_cm = np.linspace( y0_cm, y1_cm, ny + 1,)

    x_cent_cm = 0.5 * (x_edges_cm[:-1] + x_edges_cm[1:])

    y_cent_cm = 0.5 * (y_edges_cm[:-1] + y_edges_cm[1:])

    if useMicrons:
        coordinate_scale = 1.0e4
        coordinate_unit = r"$\mu$m"

        x_edges = x_edges_cm * coordinate_scale
        y_edges = y_edges_cm * coordinate_scale

        x_cent = x_cent_cm * coordinate_scale
        y_cent = y_cent_cm * coordinate_scale
    else:
        coordinate_scale = 1.0
        coordinate_unit = "cm"

        x_edges = x_edges_cm
        y_edges = y_edges_cm

        x_cent = x_cent_cm
        y_cent = y_cent_cm

    print(f"x domain: {x_edges[0]:.6g} to " f"{x_edges[-1]:.6g} {coordinate_unit}")

    print(f"y domain: {y_edges[0]:.6g} to " f"{y_edges[-1]:.6g} {coordinate_unit}")

    # ==========================================================
    # Extract requested field lineout
    # ==========================================================
    if lineout_dir == "x":
        # Fixed x, lineout versus y
        line_index = int(np.argmin(np.abs(x_cent - lineout_pos)))

        field_line = np.asarray(arr2d[line_index, :])

        line_coordinate = y_cent
        line_location = x_cent[line_index]

        line_xlabel = rf"$y$ [{coordinate_unit}]"

        line_title = (rf"{field} lineout at " rf"$x={line_location:.3f}$ {coordinate_unit}")

        # Array indexing is arr2d[x_index, y_index].
        # imshow rows therefore represent x and columns represent y.
        image_array = arr2d

        image_extent = [y_edges[0], y_edges[-1], x_edges[0], x_edges[-1],]

        image_xlabel = rf"$y$ [{coordinate_unit}]"
        image_ylabel = rf"$x$ [{coordinate_unit}]"

        marker_type = "horizontal"

    else:
        # Fixed y, lineout versus x
        line_index = int(np.argmin(np.abs(y_cent - lineout_pos)))

        field_line = np.asarray(arr2d[:, line_index])

        line_coordinate = x_cent
        line_location = y_cent[line_index]

        line_xlabel = rf"$x$ [{coordinate_unit}]"

        line_title = (rf"{field} lineout at " rf"$y={line_location:.3f}$ {coordinate_unit}")

        # Transpose so the horizontal coordinate is x.
        image_array = arr2d.T

        image_extent = [x_edges[0], x_edges[-1], y_edges[0], y_edges[-1],]

        image_xlabel = rf"$x$ [{coordinate_unit}]"
        image_ylabel = rf"$y$ [{coordinate_unit}]"

        marker_type = "horizontal"

    print(f"Requested {lineout_dir} position: " f"{lineout_pos:.6g} {coordinate_unit}")

    print(f"Nearest cell center used: " f"{line_location:.6g} {coordinate_unit}")

    # ==========================================================
    # Material-fraction lineouts
    # ==========================================================
    cham2d = np.asarray(get_field_array(cg, ftype, "cham",)[:, :, 0])

    gas2d = np.asarray(get_field_array(cg, ftype, "gas",)[:, :, 0])

    targ2d = np.asarray(get_field_array(cg, ftype, "targ",)[:, :, 0])

    expected_shape = arr2d.shape

    for name, species_array in {"cham": cham2d, "gas": gas2d, "targ": targ2d,}.items():
        if species_array.shape != expected_shape:
            raise ValueError(
                f"{name} has shape {species_array.shape}, "
                f"but {field} has shape {expected_shape}."
            )

    if lineout_dir == "x":
        chamber_line = cham2d[line_index, :]
        gas_line = gas2d[line_index, :]
        target_line = targ2d[line_index, :]
    else:
        chamber_line = cham2d[:, line_index]
        gas_line = gas2d[:, line_index]
        target_line = targ2d[:, line_index]

    # ==========================================================
    # Downsample only the displayed image
    # ==========================================================
    image_downsampled, stride = _auto_downsample_2d(image_array, max_pixels=max_plot_pixels,)

    if stride > 1:
        print(f"2D image downsampled by stride {stride}; " "lineouts use full resolution.")

    finite = image_downsampled[np.isfinite(image_downsampled)]

    finite = image_downsampled[np.isfinite(image_downsampled)]
    positive = finite[finite > 0]

    if finite.size == 0:
        raise ValueError(
            f"{field} has no finite values to plot."
        )

    # Use logarithmic scaling whenever possible
    if positive.size > 0:
        vmin = np.nanmin(positive)
        vmax = np.nanmax(finite)

        if vmax <= vmin:
            norm = Normalize(vmin=np.nanmin(finite),vmax=np.nanmax(finite),)
        else:
            norm = LogNorm(vmin=vmin,vmax=vmax,)

    else:
        norm = Normalize(vmin=np.nanmin(finite),vmax=np.nanmax(finite),)

    if field.lower() in ["velx","vely","velz"]:   # add field here for if there are negative values in field, otherwise lognorm is weird
        norm = Normalize(vmin=np.nanmin(finite),vmax=np.nanmax(finite))
    # ==========================================================
    # Figure
    # ==========================================================
    fig, (ax_2d, ax_line, ax_fraction,) = plt.subplots(3,1,figsize=(9, 11),
        gridspec_kw={"height_ratios": [3.0, 1.0, 1.0],},
        constrained_layout=True,
    )

    # ==========================================================
    # Panel 1: 2D field
    # ==========================================================
    im = ax_2d.imshow(
        image_downsampled,
        origin="lower",
        extent=image_extent,
        aspect="equal",
        interpolation="nearest",
        cmap="RdBu_r" if field.lower() in ["velx","vely","velz"] else "viridis",
        norm = norm,
        rasterized=True,
    )

    ax_2d.set_xlabel(image_xlabel)
    ax_2d.set_ylabel(image_ylabel)
    ax_2d.set_title(f"{field}, t = {sim_time_ns:.3f} ns")
    #ax_2d.set_ylim(990, 1010)
    #ax_2d.set_xlim(1980, 2010)
    if marker_type == "horizontal":
        ax_2d.axhline(line_location,linestyle="--",linewidth=1.2,color="white",)
    else:
        ax_2d.axvline(line_location,linestyle="--",linewidth=1.2,color="white",)

    colorbar = fig.colorbar(im,ax=ax_2d,fraction=0.035,pad=0.02,shrink=0.75,aspect=25,)
    colorbar.set_label(field)

    # ==========================================================
    # Panel 2: field lineout
    # ==========================================================
    ax_line.plot(line_coordinate, field_line, linewidth=2,)
    ax_line.set_ylabel(field)
    ax_line.set_title(line_title)
    ax_line.grid(True,alpha=0.3,)
    ax_line.tick_params(axis="x",labelbottom=False,)
    #ax_line.set_ylim(0,1000)
    #ax_line.set_xlim(1980, 2010)
    ax_line.ticklabel_format(axis="y",style="sci",scilimits=(0, 0),)

    # ==========================================================
    # Panel 3: material fractions
    # ==========================================================
    ax_fraction.plot(line_coordinate,chamber_line,linewidth=2,label="Chamber",)
    ax_fraction.plot(line_coordinate,gas_line,linewidth=2,label="Gas",)
    ax_fraction.plot(line_coordinate,target_line,linewidth=2,label="Target",)

    ax_fraction.set_xlabel(line_xlabel)
    ax_fraction.set_ylabel("Fraction")
    ax_fraction.set_title("Material fractions")
    ax_fraction.set_ylim(-0.05,1.05,)
    ax_fraction.grid(True,alpha=0.3,)
    ax_fraction.legend(loc="best",frameon=False,)
    #plt.show()
    # Match the horizontal ranges of the two lineout panels.
    line_min = line_coordinate[0]
    line_max = line_coordinate[-1]
    #ax_line.set_xlim(line_min,line_max,)
    #ax_fraction.set_xlim(line_min,line_max,)
    #ax_fraction.set_xlim(1980, 2010)
    # ==========================================================
    # Save
    # ==========================================================
    if savePlot:
        saveDir = Path(saveDir)

        saveDir.mkdir(parents=True,exist_ok=True,)

        unit_suffix = ("um" if useMicrons else "cm")

        output_path = saveDir / (
            f"{Path(fp).stem}_{field}_"
            f"{lineout_dir}_{lineout_pos:g}{unit_suffix}.png"
        )

        fig.savefig(output_path,dpi=150,bbox_inches="tight",)
        print(f"Saved {output_path}")

    line_data = {
        "coordinate": line_coordinate,
        "field": field_line,
        "chamber": chamber_line,
        "gas": gas_line,
        "target": target_line,
        "lineout_dir": lineout_dir,
        "requested_position": lineout_pos,
        "actual_position": line_location,
        "useMicrons": useMicrons,
        "display_stride": stride,
    }

    return (fig,(ax_2d, ax_line, ax_fraction),line_data,)

def plot_knudsen_number(
    ds, fp,
    cham={"A": 4, "Zion": 2},
    targ={"A": 20.2, "Zion": 10},
    gas={"A": 1.6, "Zion": 1.2},
    test_species="gas",
    lineout_dir="x", lineout_pos=0.0,
    useMicrons=True, savePlot=True, saveDir=Path("."),
    max_plot_pixels=1200,
    vmin=None, vmax=None, cmap="turbo"
):
    """
    Calculate and plot the effective ion-ion Knudsen number

        Kn = lambda_eff / L_T
        L_T = Ti / |grad Ti|
        1/lambda_eff = sum_j 1/lambda_ij

    Top:    2D Knudsen-number map
    Bottom: Knudsen-number lineout
    """

    # ==========================================================
    # SETTINGS
    # ==========================================================

    # ---- Figure geometry: match density plots ----
    FIGSIZE = (6.5, 8.4)
    LEFT, RIGHT = 0.19, 0.95
    BOTTOM, TOP = 0.11, 0.84
    HSPACE = 0.08
    HEIGHT_RATIOS = [3, 1.25]

    # ---- Fonts ----
    FONT_SIZE = 24
    LABEL_SIZE = 26
    TITLE_SIZE = 28
    TICK_SIZE = 24

    # ---- Lines / axes ----
    AXIS_WIDTH = 1.2
    TICK_WIDTH = 1.2
    TICK_LENGTH = 5
    LINE_WIDTH = 2.0

    # ---- Colorbar ----
    CBAR_FRACTION = 0.06
    CBAR_PAD = 0.04
    CBAR_ASPECT = 30
    CBAR_LABELPAD = 10

    plt.rcParams.update({
    # Font
    "font.family": "serif",
    "font.size": FONT_SIZE,

    # Computer Modern-style math (matches LaTeX nicely)
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

    # Lines
    "lines.linewidth": LINE_WIDTH,
    })
    # ==========================================================
    # Input setup
    # ==========================================================

    lineout_dir = lineout_dir.lower()
    test_species = test_species.lower()

    if lineout_dir not in {"x", "y"}:
        raise ValueError("lineout_dir must be either 'x' or 'y'.")

    if test_species not in {"cham", "targ", "gas"}:
        raise ValueError("test_species must be 'cham', 'targ', or 'gas'.")

    species = {"cham": cham, "targ": targ, "gas": gas}

    # ==========================================================
    # Load FLASH fields
    # ==========================================================

    cg, dims = get_covering_grid(ds)
    sim_time_ns = get_sim_time_ns(ds)

    dens2d = np.asarray(get_field_array(cg, "flash", "dens")[:, :, 0])
    tion2d = np.asarray(get_field_array(cg, "flash", "tion")[:, :, 0])
    tele2d = np.asarray(get_field_array(cg, "flash", "tele")[:, :, 0])
    cham2d = np.asarray(get_field_array(cg, "flash", "cham")[:, :, 0])
    targ2d = np.asarray(get_field_array(cg, "flash", "targ")[:, :, 0])
    gas2d  = np.asarray(get_field_array(cg, "flash", "gas")[:, :, 0])

    nx, ny = dens2d.shape
    fractions = {"cham": cham2d, "targ": targ2d, "gas": gas2d}

    # ==========================================================
    # Coordinates
    # FLASH coordinates are cm; display coordinates are mm
    # ==========================================================

    x0_cm, x1_cm = float(ds.domain_left_edge[0]), float(ds.domain_right_edge[0])
    y0_cm, y1_cm = float(ds.domain_left_edge[1]), float(ds.domain_right_edge[1])

    x_edges_cm = np.linspace(x0_cm, x1_cm, nx + 1)
    y_edges_cm = np.linspace(y0_cm, y1_cm, ny + 1)

    x_cent_cm = 0.5 * (x_edges_cm[:-1] + x_edges_cm[1:])
    y_cent_cm = 0.5 * (y_edges_cm[:-1] + y_edges_cm[1:])

    dx_cm = x_cent_cm[1] - x_cent_cm[0]
    dy_cm = y_cent_cm[1] - y_cent_cm[0]

    coordinate_scale = 10.0 if useMicrons else 1.0
    coordinate_unit = "mm" if useMicrons else "cm"

    x_edges, y_edges = x_edges_cm * coordinate_scale, y_edges_cm * coordinate_scale
    x_cent, y_cent = x_cent_cm * coordinate_scale, y_cent_cm * coordinate_scale

    # ==========================================================
    # Temperatures
    # ==========================================================

    tion_keV = tion2d * kb / (1000 * e)
    tele_keV = tele2d * kb / (1000 * e)

    # ==========================================================
    # Temperature-gradient scale length
    #
    # L_T = Ti / |grad Ti|
    # ==========================================================

    tion_smooth = gaussian_filter(tion_keV, sigma=2.0)
    dTdx, dTdy = np.gradient(tion_smooth, dx_cm, dy_cm)
    gradT = np.sqrt(dTdx**2 + dTdy**2)

    with np.errstate(divide="ignore", invalid="ignore"):
        LT_cm = np.where((gradT > 0) & np.isfinite(gradT), tion_keV / gradT, np.inf)

    LT_um = LT_cm * 1.0e4

    # ==========================================================
    # Material number densities
    #
    # n_j = rho X_j / (A_j m_p)
    # ==========================================================

    mp_g = mp * 1000.0

    rho_species = {name: dens2d * fractions[name] for name in species}
    n_cm3 = {name: rho_species[name] / (species[name]["A"] * mp_g) for name in species}
    n_1e20 = {name: n_cm3[name] / 1.0e20 for name in species}

    # ==========================================================
    # Effective ion-ion mean free path
    #
    # 1/lambda_eff = sum_j 1/lambda_ij
    # ==========================================================

    test = species[test_species]
    mi, Zi = test["A"], test["Zion"]
    lambda_ij_mm = {}

    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):

        for name in species:
            mj, Zj, nj = species[name]["A"], species[name]["Zion"], n_1e20[name]

            lambda_ij_mm[name] = mfp_ii(
                mi=mi, mj=mj, Zi=Zi, Zj=Zj,
                Ti=tion_keV, nj=nj
            )

            lambda_ij_mm[name] = np.where(nj > 0, lambda_ij_mm[name], np.inf)

        inv_lambda_mm = sum(1.0 / lambda_ij_mm[name] for name in species)
        lambda_eff_mm = np.where(inv_lambda_mm > 0, 1.0 / inv_lambda_mm, np.inf)

    lambda_eff_um = lambda_eff_mm * 1000.0

    # ==========================================================
    # Knudsen number
    # ==========================================================

    with np.errstate(divide="ignore", invalid="ignore"):
        Kn2d = lambda_eff_um / LT_um

    test_fraction = fractions[test_species]

    Kn2d = np.where(test_fraction > 0, Kn2d, np.nan)
    lambda_eff_um = np.where(test_fraction > 0, lambda_eff_um, np.nan)
    LT_um = np.where(test_fraction > 0, LT_um, np.nan)

    Kn2d = np.where(np.isfinite(Kn2d) & (Kn2d > 0), Kn2d, np.nan)

    # ==========================================================
    # Extract lineout
    # ==========================================================

    if lineout_dir == "x":

        # Hold x fixed; plot versus y
        line_index = int(np.argmin(np.abs(x_cent - lineout_pos)))

        line_coordinate = y_cent
        line_location = x_cent[line_index]

        Kn_line = Kn2d[line_index, :]
        mfp_line = lambda_eff_um[line_index, :]
        LT_line = LT_um[line_index, :]
        tion_line = tion_keV[line_index, :]
        tele_line = tele_keV[line_index, :]
        dens_line = dens2d[line_index, :]

        fraction_lines = {name: fractions[name][line_index, :] for name in species}
        density_lines = {name: n_1e20[name][line_index, :] for name in species}
        mfp_lines = {name: lambda_ij_mm[name][line_index, :] * 1000.0 for name in species}

        image_array = Kn2d
        image_extent = [y_edges[0], y_edges[-1], x_edges[0], x_edges[-1]]
        image_ylabel = rf"$x$ ({coordinate_unit})"
        line_xlabel = rf"$y$ ({coordinate_unit})"

    else:

        # Hold y fixed; plot versus x
        line_index = int(np.argmin(np.abs(y_cent - lineout_pos)))

        line_coordinate = x_cent
        line_location = y_cent[line_index]

        Kn_line = Kn2d[:, line_index]
        mfp_line = lambda_eff_um[:, line_index]
        LT_line = LT_um[:, line_index]
        tion_line = tion_keV[:, line_index]
        tele_line = tele_keV[:, line_index]
        dens_line = dens2d[:, line_index]

        fraction_lines = {name: fractions[name][:, line_index] for name in species}
        density_lines = {name: n_1e20[name][:, line_index] for name in species}
        mfp_lines = {name: lambda_ij_mm[name][:, line_index] * 1000.0 for name in species}

        image_array = Kn2d.T
        image_extent = [x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]]
        image_ylabel = rf"$y$ ({coordinate_unit})"
        line_xlabel = rf"$x$ ({coordinate_unit})"

    print(f"Requested {lineout_dir} position: {lineout_pos:.6g} {coordinate_unit}")
    print(f"Nearest cell center used: {line_location:.6g} {coordinate_unit}")

    # ==========================================================
    # Downsample image
    # ==========================================================

    image_downsampled, stride = _auto_downsample_2d(image_array, max_pixels=max_plot_pixels)

    if stride > 1:
        print(f"2D image downsampled by stride {stride}; lineout uses full resolution.")

    # ==========================================================
    # Color normalization
    # ==========================================================

    positive = image_downsampled[np.isfinite(image_downsampled) & (image_downsampled > 0)]

    if positive.size == 0:
        raise ValueError("No finite positive Knudsen numbers to plot.")

    plot_vmin = np.nanmin(positive) if vmin is None else vmin
    plot_vmax = np.nanmax(positive) if vmax is None else vmax
    norm = LogNorm(vmin=plot_vmin, vmax=plot_vmax)

    print(f"Kn color scale: {plot_vmin:.4g} -> {plot_vmax:.4g}")

    # ==========================================================
    # Figure
    # ==========================================================

    fig, (ax_2d, ax_line) = plt.subplots(
        2, 1, figsize=FIGSIZE, sharex=True,
        gridspec_kw={"height_ratios": HEIGHT_RATIOS, "hspace": HSPACE}
    )

    fig.subplots_adjust(left=LEFT, right=RIGHT, bottom=BOTTOM, top=TOP, hspace=HSPACE)

    # Apply same axis formatting as density plots
    for ax in (ax_2d, ax_line):
        for spine in ax.spines.values():
            spine.set_linewidth(AXIS_WIDTH)

        ax.tick_params(
            axis="both", which="major",
            width=TICK_WIDTH, length=TICK_LENGTH,
            labelsize=TICK_SIZE, direction="in"
        )

    # ==========================================================
    # 2D Kn map
    # ==========================================================

    im = ax_2d.imshow(
        image_downsampled, origin="lower", extent=image_extent,
        aspect="auto", interpolation="nearest",
        cmap=cmap, norm=norm, rasterized=True
    )

    ax_2d.set_ylabel(image_ylabel, fontsize=LABEL_SIZE)
    ax_2d.set_xticks([0.0, 2.5, 5.0, 7.5, 10.0])
    # Lineout location
    if lineout_dir == "x":
        ax_2d.axhline(line_location, color="white", linestyle="--", linewidth=LINE_WIDTH, alpha=0.8)
    else:
        ax_2d.axvline(line_location, color="white", linestyle="--", linewidth=LINE_WIDTH, alpha=0.8)

    # ==========================================================
    # Horizontal colorbar
    # ==========================================================

    cbar = fig.colorbar(
        im, ax=ax_2d, orientation="horizontal", location="top",
        fraction=CBAR_FRACTION, pad=CBAR_PAD, aspect=CBAR_ASPECT
    )

    cbar.ax.tick_params(width=TICK_WIDTH, length=TICK_LENGTH, labelsize=TICK_SIZE, pad=3)
    cbar.outline.set_linewidth(AXIS_WIDTH)
    cbar.set_label(r"$K_{n}$", fontsize=LABEL_SIZE, labelpad=CBAR_LABELPAD)

    # ==========================================================
    # Kn lineout
    # ==========================================================

    ax_line.plot(line_coordinate, Kn_line, linewidth=LINE_WIDTH)
    ax_line.axhline(1.0, color="black", linestyle="--", linewidth=1.2, alpha=0.7)

    ax_line.set_yscale("log")
    ax_line.set_xlabel(line_xlabel, fontsize=LABEL_SIZE)
    ax_line.set_ylabel(r"$K_{n}$", fontsize=LABEL_SIZE)
    ax_line.set_xticks([0.0, 2.5, 5.0, 7.5, 10.0])
    ax_line.set_yticks([1e-7, 1e-4, 1e-1])
    ax_line.grid(True, alpha=0.25)

    # ==========================================================
    # Save
    # ==========================================================

    output_path = None

    if savePlot:
        saveDir = Path(saveDir)
        saveDir.mkdir(parents=True, exist_ok=True)

        unit_suffix = "mm" if useMicrons else "cm"

        output_path = saveDir / (
            f"{Path(fp).stem}_Kn_{test_species}_"
            f"{lineout_dir}_{lineout_pos:g}{unit_suffix}.pdf"
        )

        fig.savefig(output_path)
        print(f"Saved {output_path}")

    plt.show()

    # ==========================================================
    # Return data
    # ==========================================================

    line_data = {
        "coordinate": line_coordinate,
        "Kn": Kn_line,
        "lambda_eff_um": mfp_line,
        "LT_um": LT_line,
        "Ti_keV": tion_line,
        "Te_keV": tele_line,
        "density_gcm3": dens_line,
        "fractions": fraction_lines,
        "n_1e20": density_lines,
        "lambda_ij_um": mfp_lines,
        "test_species": test_species,
        "actual_position": line_location,
        "display_stride": stride,
        "output_path": output_path,
    }

    return fig, (ax_2d, ax_line), line_data

def plot_2d_ne_with_material_fractions(
    ds,
    fp,
    ftype,
    lineout_dir="x",
    lineout_pos=0.0,
    useMicrons=True,
    savePlot=True,
    saveDir=Path("."),
    max_plot_pixels=1200,
):
    """
    Plot 2D electron number density, electron/ion density lineouts,
    and material fractions.

    Electron density:
        ne = rho * ye / mp

    Ion density:
        ni = rho / (abar * mp)

    Parameters
    ----------
    ds
        Loaded yt FLASH dataset.

    fp
        Path to the FLASH plotfile.

    ftype
        FLASH field type passed to get_field_array().

    lineout_dir : {"x", "y"}
        Coordinate held fixed.

        lineout_dir="x":
            Hold x fixed at lineout_pos and plot versus y.

        lineout_dir="y":
            Hold y fixed at lineout_pos and plot versus x.

    lineout_pos : float
        In microns when useMicrons=True.
        In cm when useMicrons=False.

    useMicrons : bool
        Convert plotted spatial coordinates from cm to microns.

    savePlot : bool
        Save the figure.

    saveDir
        Directory where figure is saved.

    Returns
    -------
    fig, axes, line_data
    """

    lineout_dir = lineout_dir.lower()

    if lineout_dir not in {"x", "y"}:
        raise ValueError(
            "lineout_dir must be either 'x' or 'y'."
        )

    # ==========================================================
    # Covering grid
    # ==========================================================
    cg, dims = get_covering_grid(ds)
    sim_time_ns = get_sim_time_ns(ds)

    # ----------------------------------------------------------
    # FLASH quantities
    # ----------------------------------------------------------
    rho = get_field_array(cg, ftype, "dens")
    ye = get_field_array(cg, ftype, "ye")
    abar = get_field_array(cg, ftype, "gas")

    # ----------------------------------------------------------
    # Number densities [cm^-3]
    #
    # rho  : g/cm^3
    # mp   : g
    # ye   : electrons / baryon
    # abar : mean atomic mass in proton-mass units
    # ----------------------------------------------------------
    mp = constants.proton_mass # mass of H(kg)
    mp = mp*1000
    ne = (rho / mp) * ye

    ni = rho / (abar * mp)

    ne2d = np.asarray(ne[:, :, 0])
    ni2d = np.asarray(ni[:, :, 0])

    if ne2d.ndim != 2:
        raise ValueError(
            f"Expected a 2D electron density array, "
            f"got shape {ne2d.shape}."
        )

    if ni2d.shape != ne2d.shape:
        raise ValueError(
            f"ni has shape {ni2d.shape}, "
            f"but ne has shape {ne2d.shape}."
        )

    nx, ny = ne2d.shape

    print("ne shape:", ne2d.shape)

    # ==========================================================
    # Coordinate arrays
    # ==========================================================
    x0_cm = float(ds.domain_left_edge[0])
    x1_cm = float(ds.domain_right_edge[0])

    y0_cm = float(ds.domain_left_edge[1])
    y1_cm = float(ds.domain_right_edge[1])

    x_edges_cm = np.linspace(
        x0_cm, x1_cm, nx + 1
    )

    y_edges_cm = np.linspace(
        y0_cm, y1_cm, ny + 1
    )

    x_cent_cm = 0.5 * (
        x_edges_cm[:-1] + x_edges_cm[1:]
    )

    y_cent_cm = 0.5 * (
        y_edges_cm[:-1] + y_edges_cm[1:]
    )

    if useMicrons:

        coordinate_scale = 1.0e4
        coordinate_unit = r"$\mu$m"

    else:

        coordinate_scale = 1.0
        coordinate_unit = "cm"

    x_edges = x_edges_cm * coordinate_scale
    y_edges = y_edges_cm * coordinate_scale

    x_cent = x_cent_cm * coordinate_scale
    y_cent = y_cent_cm * coordinate_scale

    # ==========================================================
    # Extract ne and ni lineouts
    # ==========================================================
    if lineout_dir == "x":

        # Hold x fixed; lineout versus y
        line_index = int(
            np.argmin(np.abs(x_cent - lineout_pos))
        )

        ne_line = ne2d[line_index, :]
        ni_line = ni2d[line_index, :]

        line_coordinate = y_cent
        line_location = x_cent[line_index]

        line_xlabel = rf"$y$ [{coordinate_unit}]"

        line_title = (
            rf"Density lineout at "
            rf"$x={line_location:.3f}$ {coordinate_unit}"
        )

        image_array = ne2d

        image_extent = [
            y_edges[0],
            y_edges[-1],
            x_edges[0],
            x_edges[-1],
        ]

        image_xlabel = rf"$y$ [{coordinate_unit}]"
        image_ylabel = rf"$x$ [{coordinate_unit}]"

    else:

        # Hold y fixed; lineout versus x
        line_index = int(
            np.argmin(np.abs(y_cent - lineout_pos))
        )

        ne_line = ne2d[:, line_index]
        ni_line = ni2d[:, line_index]

        line_coordinate = x_cent
        line_location = y_cent[line_index]

        line_xlabel = rf"$x$ [{coordinate_unit}]"

        line_title = (
            rf"Density lineout at "
            rf"$y={line_location:.3f}$ {coordinate_unit}"
        )

        image_array = ne2d.T

        image_extent = [
            x_edges[0],
            x_edges[-1],
            y_edges[0],
            y_edges[-1],
        ]

        image_xlabel = rf"$x$ [{coordinate_unit}]"
        image_ylabel = rf"$y$ [{coordinate_unit}]"

    print(
        f"Requested {lineout_dir} position: "
        f"{lineout_pos:.6g} {coordinate_unit}"
    )

    print(
        f"Nearest cell center used: "
        f"{line_location:.6g} {coordinate_unit}"
    )

    # ==========================================================
    # Material fractions
    # ==========================================================
    cham2d = np.asarray(
        get_field_array(cg, ftype, "cham")[:, :, 0]
    )

    gas2d = np.asarray(
        get_field_array(cg, ftype, "gas")[:, :, 0]
    )

    targ2d = np.asarray(
        get_field_array(cg, ftype, "targ")[:, :, 0]
    )

    if lineout_dir == "x":

        chamber_line = cham2d[line_index, :]
        gas_line = gas2d[line_index, :]
        target_line = targ2d[line_index, :]

    else:

        chamber_line = cham2d[:, line_index]
        gas_line = gas2d[:, line_index]
        target_line = targ2d[:, line_index]

    # ==========================================================
    # Downsample ONLY the 2D displayed image
    # ==========================================================
    image_downsampled, stride = _auto_downsample_2d(
        image_array,
        max_pixels=max_plot_pixels,
    )

    if stride > 1:
        print(
            f"2D image downsampled by stride {stride}; "
            "lineouts use full resolution."
        )

    finite = image_downsampled[
        np.isfinite(image_downsampled)
    ]

    positive = finite[finite > 0]

    if positive.size == 0:
        raise ValueError(
            "Electron density contains no positive finite values."
        )

    norm = LogNorm(
        vmin=np.nanmin(positive),
        vmax=np.nanmax(positive),
    )

    # ==========================================================
    # Figure
    # ==========================================================
    fig, (
        ax_2d,
        ax_line,
        ax_fraction,
    ) = plt.subplots(
        3,
        1,
        figsize=(9, 11),
        gridspec_kw={
            "height_ratios": [3.0, 1.0, 1.0]
        },
        constrained_layout=True,
    )

    # ==========================================================
    # Panel 1: 2D electron density
    # ==========================================================
    im = ax_2d.imshow(
        image_downsampled,
        origin="lower",
        extent=image_extent,
        aspect="equal",
        interpolation="nearest",
        cmap="viridis",
        norm=norm,
        rasterized=True,
    )

    ax_2d.set_xlabel(image_xlabel)
    ax_2d.set_ylabel(image_ylabel)

    ax_2d.set_title(
        rf"$n_e$, t = {sim_time_ns:.3f} ns"
    )

    # Your lineout convention always corresponds to
    # a horizontal line in the displayed orientation.
    ax_2d.axhline(
        line_location,
        linestyle="--",
        linewidth=1.2,
        color="white",
    )

    colorbar = fig.colorbar(
        im,
        ax=ax_2d,
        fraction=0.035,
        pad=0.02,
        shrink=0.75,
        aspect=25,
    )

    colorbar.set_label(
        r"$n_e$ [cm$^{-3}$]"
    )

    # ==========================================================
    # Panel 2: ne and ni lineout
    # ==========================================================
    ax_line.plot(
        line_coordinate,
        ne_line,
        linewidth=2,
        label=r"$n_e$",
    )

    # ax_line.plot(
    #     line_coordinate,
    #     ni_line,
    #     linewidth=2,
    #     label=r"$n_i$",
    # )

    ax_line.set_ylabel(
        r"Number density [cm$^{-3}$]"
    )

    ax_line.set_title(line_title)

    ax_line.grid(
        True,
        alpha=0.3,
    )

    ax_line.legend(
        loc="best",
        frameon=False,
    )

    ax_line.tick_params(
        axis="x",
        labelbottom=False,
    )

    ax_line.ticklabel_format(
        axis="y",
        style="sci",
        scilimits=(0, 0),
    )

    # ==========================================================
    # Panel 3: material fractions
    # ==========================================================
    ax_fraction.plot(
        line_coordinate,
        chamber_line,
        linewidth=2,
        label="Chamber",
    )

    ax_fraction.plot(
        line_coordinate,
        gas_line,
        linewidth=2,
        label="Gas",
    )

    ax_fraction.plot(
        line_coordinate,
        target_line,
        linewidth=2,
        label="Target",
    )

    ax_fraction.set_xlabel(line_xlabel)
    ax_fraction.set_ylabel("Fraction")
    ax_fraction.set_title("Material fractions")

    ax_fraction.set_ylim(
        -0.05,
        1.05,
    )

    ax_fraction.grid(
        True,
        alpha=0.3,
    )

    ax_fraction.legend(
        loc="best",
        frameon=False,
    )

    # Match horizontal ranges
    line_min = line_coordinate[0]
    line_max = line_coordinate[-1]

    ax_line.set_xlim(
        line_min,
        line_max,
    )

    ax_fraction.set_xlim(
        line_min,
        line_max,
    )
    plt.show()
    # ==========================================================
    # Save
    # ==========================================================
    if savePlot:

        saveDir = Path(saveDir)

        saveDir.mkdir(
            parents=True,
            exist_ok=True,
        )

        unit_suffix = (
            "um" if useMicrons else "cm"
        )

        output_path = saveDir / (
            f"{Path(fp).stem}_ne_"
            f"{lineout_dir}_"
            f"{lineout_pos:g}{unit_suffix}.png"
        )

        fig.savefig(
            output_path,
            dpi=150,
            bbox_inches="tight",
        )

        print(f"Saved {output_path}")

    # ==========================================================
    # Returned lineout data
    # ==========================================================
    line_data = {
        "coordinate": line_coordinate,
        "ne": np.asarray(ne_line),
        "ni": np.asarray(ni_line),
        "chamber": chamber_line,
        "gas": gas_line,
        "target": target_line,
        "lineout_dir": lineout_dir,
        "requested_position": lineout_pos,
        "actual_position": line_location,
        "useMicrons": useMicrons,
        "display_stride": stride,
    }

    return (
        fig,
        (ax_2d, ax_line, ax_fraction),
        line_data,
    )


def plot_2d_fields(
    ds,
    fp,
    field,
    ftype="flash",
    lineout_dir="x",
    lineout_pos=0.0,
    useMicrons=True,
    savePlots=True,
    saveDir=Path("."),
    max_plot_pixels=1200,
):
    """
    Plot a general 2D FLASH field with a full-resolution lineout.

    Parameters
    ----------
    ds
        Loaded yt FLASH dataset.

    fp
        Path to FLASH plotfile.

    ftype
        FLASH field type, usually "flash".

    field
        FLASH field name, e.g. "dens", "pres", "tele", "velx".

    lineout_dir : {"x", "y"}
        Coordinate held fixed.

        lineout_dir="x":
            Hold x fixed at lineout_pos and plot versus y.

        lineout_dir="y":
            Hold y fixed at lineout_pos and plot versus x.

    lineout_pos : float
        Requested fixed-coordinate position.

        Microns if useMicrons=True.
        cm if useMicrons=False.

    useMicrons : bool
        Convert coordinates from cm to microns.

    savePlot : bool
        Save figure if True.

    saveDir
        Output directory.

    max_plot_pixels : int
        Approximate maximum number of pixels along longest
        displayed image dimension.

        Only the displayed 2D image is downsampled.
        The lineout remains full resolution.

    Returns
    -------
    fig, axes, line_data

        axes = (ax_2d, ax_line)

        line_data contains:
            coordinate
            field
            requested_position
            actual_position
            lineout_dir
            display_stride
    """

    # ==========================================================
    # Check lineout direction
    # ==========================================================
    lineout_dir = lineout_dir.lower()

    if lineout_dir not in {"x", "y"}:
        raise ValueError(
            "lineout_dir must be either 'x' or 'y'."
        )

    # ==========================================================
    # Load covering grid and requested field
    # ==========================================================
    cg, dims = get_covering_grid(ds)

    sim_time_ns = get_sim_time_ns(ds)

    # IMPORTANT:
    # Do not automatically multiply by kb/e here.
    # We want this function to work for ANY FLASH field.
    arr = get_field_array(
        cg,
        ftype,
        field,
    )

    arr2d = np.asarray(
        arr[:, :, 0]
    )

    if arr2d.ndim != 2:
        raise ValueError(
            f"Expected a 2D field array, got shape {arr2d.shape}."
        )

    nx, ny = arr2d.shape

    print(
        f"{field} shape: {arr2d.shape}"
    )

    # ==========================================================
    # Coordinate arrays
    # ==========================================================
    x0_cm = float(ds.domain_left_edge[0])
    x1_cm = float(ds.domain_right_edge[0])

    y0_cm = float(ds.domain_left_edge[1])
    y1_cm = float(ds.domain_right_edge[1])

    # Cell edges
    x_edges_cm = np.linspace(
        x0_cm,
        x1_cm,
        nx + 1,
    )

    y_edges_cm = np.linspace(
        y0_cm,
        y1_cm,
        ny + 1,
    )

    # Cell centers
    x_cent_cm = 0.5 * (
        x_edges_cm[:-1]
        + x_edges_cm[1:]
    )

    y_cent_cm = 0.5 * (
        y_edges_cm[:-1]
        + y_edges_cm[1:]
    )

    # ==========================================================
    # Coordinate units
    # ==========================================================
    if useMicrons:

        coordinate_scale = 1.0e4
        coordinate_unit = r"$\mu$m"

    else:

        coordinate_scale = 1.0
        coordinate_unit = "cm"

    x_edges = x_edges_cm * coordinate_scale
    y_edges = y_edges_cm * coordinate_scale

    x_cent = x_cent_cm * coordinate_scale
    y_cent = y_cent_cm * coordinate_scale

    print(
        f"x domain: {x_edges[0]:.6g} to "
        f"{x_edges[-1]:.6g} {coordinate_unit}"
    )

    print(
        f"y domain: {y_edges[0]:.6g} to "
        f"{y_edges[-1]:.6g} {coordinate_unit}"
    )

    # ==========================================================
    # Extract requested lineout
    # ==========================================================
    if lineout_dir == "x":

        # --------------------------------------------
        # Hold x fixed; lineout runs along y
        # --------------------------------------------
        line_index = int(
            np.argmin(
                np.abs(
                    x_cent - lineout_pos
                )
            )
        )

        field_line = np.asarray(
            arr2d[line_index, :]
        )

        line_coordinate = y_cent

        line_location = x_cent[
            line_index
        ]

        line_xlabel = (
            rf"$y$ [{coordinate_unit}]"
        )

        line_title = (
            rf"{field} lineout at "
            rf"$x={line_location:.3f}$ "
            rf"{coordinate_unit}"
        )

        # arr2d indexing:
        #
        # arr2d[x_index, y_index]
        #
        # imshow rows therefore correspond to x
        # and columns correspond to y.
        image_array = arr2d

        image_extent = [
            y_edges[0],
            y_edges[-1],
            x_edges[0],
            x_edges[-1],
        ]

        image_xlabel = (
            rf"$y$ [{coordinate_unit}]"
        )

        image_ylabel = (
            rf"$x$ [{coordinate_unit}]"
        )

        marker_type = "horizontal"

    else:

        # --------------------------------------------
        # Hold y fixed; lineout runs along x
        # --------------------------------------------
        line_index = int(
            np.argmin(
                np.abs(
                    y_cent - lineout_pos
                )
            )
        )

        field_line = np.asarray(
            arr2d[:, line_index]
        )

        line_coordinate = x_cent

        line_location = y_cent[
            line_index
        ]

        line_xlabel = (
            rf"$x$ [{coordinate_unit}]"
        )

        line_title = (
            rf"{field} lineout at "
            rf"$y={line_location:.3f}$ "
            rf"{coordinate_unit}"
        )

        # Transpose so horizontal plot coordinate is x
        image_array = arr2d.T

        image_extent = [
            x_edges[0],
            x_edges[-1],
            y_edges[0],
            y_edges[-1],
        ]

        image_xlabel = (
            rf"$x$ [{coordinate_unit}]"
        )

        image_ylabel = (
            rf"$y$ [{coordinate_unit}]"
        )

        marker_type = "horizontal"

    print(
        f"Requested {lineout_dir} position: "
        f"{lineout_pos:.6g} {coordinate_unit}"
    )

    print(
        f"Nearest cell center used: "
        f"{line_location:.6g} {coordinate_unit}"
    )

    # ==========================================================
    # Downsample displayed image only
    # ==========================================================
    image_downsampled, stride = (
        _auto_downsample_2d(
            image_array,
            max_pixels=max_plot_pixels,
        )
    )

    if stride > 1:

        print(
            f"2D image downsampled by stride {stride}; "
            "lineout uses full resolution."
        )

    # ==========================================================
    # Determine normalization
    # ==========================================================
    finite = image_downsampled[
        np.isfinite(
            image_downsampled
        )
    ]

    if finite.size == 0:
        raise ValueError(
            f"{field} has no finite values to plot."
        )

    # Signed quantities should NOT use LogNorm.
    signed_fields = {
        "velx",
        "vely",
        "velz",
    }

    if field.lower() in signed_fields:

        norm = Normalize(
            vmin=np.nanmin(finite),
            vmax=np.nanmax(finite),
        )

    else:

        positive = finite[
            finite > 0
        ]

        if positive.size > 0:

            vmin = np.nanmin(
                positive
            )

            vmax = np.nanmax(
                finite
            )

            if vmax > vmin:

                norm = LogNorm(
                    vmin=vmin,
                    vmax=vmax,
                )

            else:

                norm = Normalize(
                    vmin=np.nanmin(finite),
                    vmax=np.nanmax(finite),
                )

        else:

            norm = Normalize(
                vmin=np.nanmin(finite),
                vmax=np.nanmax(finite),
            )

    # ==========================================================
    # Figure
    # ==========================================================
    fig, (
        ax_2d,
        ax_line,
    ) = plt.subplots(
        2,
        1,
        figsize=(9, 8),
        gridspec_kw={
            "height_ratios": [
                3.0,
                1.0,
            ]
        },
        constrained_layout=True,
    )

    # ==========================================================
    # Panel 1: 2D field
    # ==========================================================
    im = ax_2d.imshow(
        image_downsampled,
        origin="lower",
        extent=image_extent,
        aspect="equal",
        interpolation="nearest",
        cmap="viridis",
        norm=norm,
        rasterized=True,
    )

    ax_2d.set_xlabel(
        image_xlabel
    )

    ax_2d.set_ylabel(
        image_ylabel
    )

    ax_2d.set_title(
        f"{field}, "
        f"t = {sim_time_ns:.3f} ns"
    )

    # Draw chosen lineout position
    if marker_type == "horizontal":

        ax_2d.axhline(
            line_location,
            linestyle="--",
            linewidth=1.2,
            color="white",
        )

    else:

        ax_2d.axvline(
            line_location,
            linestyle="--",
            linewidth=1.2,
            color="white",
        )

    # Colorbar
    colorbar = fig.colorbar(
        im,
        ax=ax_2d,
        fraction=0.035,
        pad=0.02,
        shrink=0.75,
        aspect=25,
    )

    colorbar.set_label(
        field
    )

    # ==========================================================
    # Panel 2: lineout
    # ==========================================================
    ax_line.plot(
        line_coordinate,
        field_line,
        linewidth=2,
    )

    ax_line.set_xlabel(
        line_xlabel
    )

    ax_line.set_ylabel(
        field
    )

    ax_line.set_title(
        line_title
    )

    ax_line.grid(
        True,
        alpha=0.3,
    )

    ax_line.ticklabel_format(
        axis="y",
        style="sci",
        scilimits=(0, 0),
    )

    ax_line.set_xlim(
        line_coordinate[0],
        line_coordinate[-1],
    )

    # ==========================================================
    # Save
    # ==========================================================
    if savePlots:

        saveDir = Path(
            saveDir
        )

        saveDir.mkdir(
            parents=True,
            exist_ok=True,
        )

        unit_suffix = (
            "um"
            if useMicrons
            else "cm"
        )

        output_path = (
            saveDir
            / (
                f"{Path(fp).stem}_"
                f"{field}_"
                f"{lineout_dir}_"
                f"{lineout_pos:g}"
                f"{unit_suffix}.png"
            )
        )

        fig.savefig(
            output_path,
            dpi=150,
            bbox_inches="tight",
        )

        print(
            f"Saved {output_path}"
        )

    # ==========================================================
    # Return lineout data
    # ==========================================================
    line_data = {
        "coordinate": line_coordinate,
        "field": field_line,
        "lineout_dir": lineout_dir,
        "requested_position": lineout_pos,
        "actual_position": line_location,
        "useMicrons": useMicrons,
        "display_stride": stride,
    }

    return (
        fig,
        (
            ax_2d,
            ax_line,
        ),
        line_data,
    )


def read_flash_rays(filename):
    """
    Reads FLASH RayData from a plot/checkpoint file.

    Assumes columns:
        0 = ray tag
        1 = x
        2 = y
        3 = z
        4 = power
    """
    with h5py.File(filename, "r") as f:
        if "RayData" not in f:
            raise ValueError(f"No RayData found in {filename}")

        ray_data = f["RayData"][:]

    return ray_data


def plot_rays(
    ax,
    ray_data,
    xscale=1.0,
    yscale=1.0,
    color="w",
    lw=1.0,
    cylindrical=False,
):
    """
    Overlay FLASH rays on an existing matplotlib axis.

    Parameters
    ----------
    ax
        Matplotlib axis to plot on.

    ray_data
        FLASH ray data array.

        Assumed columns:
            0 : ray tag / ID
            1 : x [cm]
            2 : z [cm]
            3 : y [cm]
            4 : ray power / weight

    xscale, yscale
        Unit conversion factors.
        Example: FLASH cm -> micron means scale = 1e4.

    cylindrical : bool, default=False
        If False:
            Plot the original coordinates:
                horizontal = column 1
                vertical   = column 2

        If True:
            Project the 3D ray coordinates into cylindrical r-z:
                horizontal = z
                vertical   = r = sqrt(x^2 + y^2)
    """

    tags = ray_data[:, 0]
    unique_tags = np.unique(tags)

    for tag in unique_tags:

        this_ray = ray_data[tags == tag]

        if cylindrical:

            # Native FLASH 3D coordinates
            ray_x = this_ray[:, 1]
            ray_z = this_ray[:, 2]
            ray_y = this_ray[:, 3]

            # Cylindrical radius
            ray_r = np.sqrt(ray_x**2 + ray_y**2)

            # Cylindrical r-z plotting convention:
            # horizontal = z
            # vertical   = r
            x = ray_z * xscale
            y = ray_r * yscale

        else:

            # Preserve original behavior
            x = this_ray[:, 1] * xscale
            y = this_ray[:, 2] * yscale

        ax.plot(
            x,
            y,
            color=color,
            lw=lw,
            alpha=0.8,
        )


def make_ray_diagnostic_figure(
    ds,
    fp,
    ray_data,
    geometry="cartesian",   # "cartesian" or "cylindrical"
    useMicrons=True,
    savePlots=False,
    saveDir=Path("."),
):
    """
    Ray diagnostic figure:
      - histogram of stored initial ray powers
      - title/text box contains:
          number of stored rays
          sampled total initial power
          min/max ray power
          total deposited energy

    Assumes:
      ray_data[:, 0] = ray tag
      ray_data[:, 4] = ray power in erg/s
      depo has units erg/g
      dens has units g/cm^3
    """

    cg, dims = get_covering_grid(ds)
    sim_time_ns = get_sim_time_ns(ds)

    dens = get_field_array(cg, "flash", "dens")[:, :, 0]
    depo = get_field_array(cg, "flash", "depo")[:, :, 0]

    nx, ny = dens.shape

    x0, x1 = float(ds.domain_left_edge[0]), float(ds.domain_right_edge[0])
    y0, y1 = float(ds.domain_left_edge[1]), float(ds.domain_right_edge[1])

    x_edges = np.linspace(x0, x1, nx + 1)
    y_edges = np.linspace(y0, y1, ny + 1)

    dx = np.diff(x_edges)
    dy = np.diff(y_edges)

    # ------------------------------------------------------------
    # Volume / area weighting
    # ------------------------------------------------------------
    if geometry.lower() in ["cartesian", "cart"]:
        # 2D Cartesian: this is energy per unit depth.
        # Units: erg/cm if depo is erg/g and dens is g/cm^3.
        dV = dx[:, None] * dy[None, :]
        Edep_erg = np.nansum(dens * depo * dV)
        Edep_label = r"$E_{\rm dep}$ per unit depth"

    elif geometry.lower() in ["cylindrical", "cyl", "rz"]:
        # 2D cylindrical R-Z:
        # x is r, y is z
        r_cent = 0.5 * (x_edges[:-1] + x_edges[1:])

        dV = 2.0 * np.pi * r_cent[:, None] * dx[:, None] * dy[None, :]
        Edep_erg = np.nansum(dens * depo * dV)
        Edep_label = r"$E_{\rm dep}$"

    else:
        raise ValueError("geometry must be 'cartesian' or 'cylindrical'.")

    Edep_J = Edep_erg * 1e-7

    # ------------------------------------------------------------
    # Ray powers
    # ------------------------------------------------------------
    tags = ray_data[:, 0]
    unique_tags = np.unique(tags)

    ray_initial_power_W = []

    for tag in unique_tags:
        r = ray_data[tags == tag]
        ray_initial_power_W.append(r[0, 4] * 1e-7)  # erg/s -> W

    ray_initial_power_W = np.array(ray_initial_power_W)

    n_rays = len(unique_tags)
    sampled_power_TW = np.nansum(ray_initial_power_W) / 1e12
    pmin = np.nanmin(ray_initial_power_W)
    pmax = np.nanmax(ray_initial_power_W)

    # ------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(7, 5))

    ax.hist(ray_initial_power_W, bins=40)

    ax.set_xlabel("Initial stored ray power [W]")
    ax.set_ylabel("Number of stored rays")

    if geometry.lower() in ["cartesian", "cart"]:
        title_geom = "2D Cartesian"
    else:
        title_geom = "2D Cylindrical R-Z"

    ax.set_title(
        f"{title_geom} ray diagnostic, t = {sim_time_ns:.3f} ns\n"
        f"{Edep_label} = {Edep_J:.3e} J"
    )

    info_text = (
        f"stored rays = {n_rays}\n"
        f"sampled initial power = {sampled_power_TW:.3f} TW\n"
        f"ray power min = {pmin:.3e} W\n"
        f"ray power max = {pmax:.3e} W"
    )

    ax.text(
        0.98,
        0.95,
        info_text,
        transform=ax.transAxes,
        ha="right",
        va="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
    )

    ax.grid(True, alpha=0.3)

    plt.tight_layout()

    if savePlots:
        saveDir = Path(saveDir)
        saveDir.mkdir(parents=True, exist_ok=True)

        out = saveDir / f"{Path(fp).stem}_ray_diagnostic_{geometry}.png"
        plt.savefig(out, dpi=200, bbox_inches="tight")
        print(f"Saved {out}")

    plt.show()

    return {
        "fig": fig,
        "ax": ax,
        "Edep_erg": Edep_erg,
        "Edep_J": Edep_J,
        "n_stored_rays": n_rays,
        "sampled_power_TW": sampled_power_TW,
        "ray_power_W": ray_initial_power_W,
        "ray_power_min_W": pmin,
        "ray_power_max_W": pmax,
    }


def plot_density_lineout(
    ds,
    fp,
    ftype,
    savePlots=False,
    saveDir=Path("."),
):
    """
    Plot a 1D FLASH density lineout with material-interface markers.

    Hardcoded lineout choices:
      line_axis = "x" -> hold x fixed and plot versus y
      line_axis = "y" -> hold y fixed and plot versus x

    Coordinates are always plotted in microns.
    """

    # -----------------------
    # Hardcoded options
    # -----------------------
    line_axis = "x"       # "x" or "y"
    line_value_um = 0
    save_dpi = 150

    # -----------------------
    # Load FLASH data
    # -----------------------
    cg, _ = get_covering_grid(ds)
    sim_time_ns = get_sim_time_ns(ds)

    dens2d = np.asarray(
        get_field_array(cg, ftype, "dens")[:, :, 0]
    )

    cham2d = np.asarray(
        get_field_array(cg, ftype, "cham")[:, :, 0]
    )
    gas2d = np.asarray(
        get_field_array(cg, ftype, "gas")[:, :, 0]
    )
    targ2d = np.asarray(
        get_field_array(cg, ftype, "targ")[:, :, 0]
    )

    # -----------------------
    # Coordinates in microns
    # -----------------------
    x0 = float(ds.domain_left_edge[0])
    x1 = float(ds.domain_right_edge[0])
    y0 = float(ds.domain_left_edge[1])
    y1 = float(ds.domain_right_edge[1])

    nx, ny = dens2d.shape

    x_edges_um = np.linspace(x0, x1, nx + 1) * 1e4
    y_edges_um = np.linspace(y0, y1, ny + 1) * 1e4

    x_cent_um = 0.5 * (
        x_edges_um[:-1] + x_edges_um[1:]
    )
    y_cent_um = 0.5 * (
        y_edges_um[:-1] + y_edges_um[1:]
    )

    # -----------------------
    # Extract lineouts
    # -----------------------
    if line_axis.lower() == "x":
        index = np.argmin(
            np.abs(x_cent_um - line_value_um)
        )

        density_line = dens2d[index, :]
        cham_line = cham2d[index, :]
        gas_line = gas2d[index, :]
        targ_line = targ2d[index, :]

        line_coord = y_cent_um
        line_location = x_cent_um[index]

        coordinate_label = r"$y$ [$\mu$m]"
        title = (
            rf"Density lineout at "
            rf"$x={line_location:.3f}\ \mu$m, "
            rf"$t={sim_time_ns:.3f}$ ns"
        )

    elif line_axis.lower() == "y":
        index = np.argmin(
            np.abs(y_cent_um - line_value_um)
        )

        density_line = dens2d[:, index]
        cham_line = cham2d[:, index]
        gas_line = gas2d[:, index]
        targ_line = targ2d[:, index]

        line_coord = x_cent_um
        line_location = y_cent_um[index]

        coordinate_label = r"$x$ [$\mu$m]"
        title = (
            rf"Density lineout at "
            rf"$y={line_location:.3f}\ \mu$m, "
            rf"$t={sim_time_ns:.3f}$ ns"
        )

    else:
        raise ValueError(
            "line_axis must be either 'x' or 'y'."
        )

    # -----------------------
    # Find material interfaces
    # -----------------------
    fractions = np.vstack(
        [cham_line, gas_line, targ_line]
    )

    # 0 = chamber, 1 = gas, 2 = target
    dominant_species = np.argmax(
        fractions,
        axis=0,
    )

    transition_indices = np.where(
        np.diff(dominant_species) != 0
    )[0]

    species_names = {
        0: "Chamber",
        1: "Gas",
        2: "Target",
    }

    interfaces = []

    for i in transition_indices:
        position = 0.5 * (
            line_coord[i] + line_coord[i + 1]
        )

        left_name = species_names[
            dominant_species[i]
        ]
        right_name = species_names[
            dominant_species[i + 1]
        ]

        interfaces.append(
            (position, f"{left_name}/{right_name}")
        )

    print(
        f"Density min/max along lineout: "
        f"{np.nanmin(density_line):.6e}, "
        f"{np.nanmax(density_line):.6e}"
    )

    print("Material interfaces:")
    for position, label in interfaces:
        print(f"  {label}: {position:.4f} microns")

    # -----------------------
    # Plot
    # -----------------------
    fig, ax = plt.subplots(
        figsize=(8, 4.5),
        constrained_layout=True,
    )

    ax.plot(
        line_coord,
        density_line,
        linewidth=2,
        label="Density",
    )

    for position, label in interfaces:
        ax.axvline(
            position,
            linestyle="--",
            linewidth=1.5,
            label=label,
        )

    ax.set_xlabel(coordinate_label)
    ax.set_ylabel(r"$\rho$ [g/cm$^3$]")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)

    ax.ticklabel_format(
        axis="y",
        style="sci",
        scilimits=(0, 0),
    )

    if interfaces:
        # Remove duplicate labels from the legend
        handles, labels = ax.get_legend_handles_labels()
        unique = dict(zip(labels, handles))

        ax.legend(
            unique.values(),
            unique.keys(),
            loc="best",
        )
    plt.show()
    # -----------------------
    # Save
    # -----------------------
    if savePlots:
        saveDir = Path(saveDir)
        saveDir.mkdir(
            parents=True,
            exist_ok=True,
        )

        output_path = saveDir / (
            f"{Path(fp).stem}_density_lineout_"
            f"{line_axis}{line_value_um:.1f}um.png"
        )

        fig.savefig(
            output_path,
            dpi=save_dpi,
            bbox_inches="tight",
        )

        print(f"Saved {output_path}")

    return (
        fig,
        ax,
        line_coord,
        density_line,
        interfaces,
    )