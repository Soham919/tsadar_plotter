import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from pathlib import Path
import yt
import glob
import os
import re

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


def plot_1d_profiles(ds, fields, sim_time_ns, useMicrons, savePlots, saveDir, fp, rays = False):
    """
    Plot 1D profiles from FLASH output.
    """
    cg, dims = get_covering_grid(ds)
    sim_time_ns = get_sim_time_ns(ds)
    x0 = float(ds.domain_left_edge[0])
    x1 = float(ds.domain_right_edge[0])

    # Get nx from any valid field
    example_field = None
    for f in fields:
        if ("flash", f) in ds.field_list or ("flash", f) in ds.derived_field_list:
            example_field = f
            break

    if example_field is None:
        raise ValueError("None of the requested 1D fields were found in the dataset.")

    example_arr = squeeze_to_1d(get_field_array(cg, example_field))
    nx = example_arr.size

    # Cell centers, not edges
    x_cm = np.linspace(x0, x1, nx)
    x_plot, xlabel = coord_to_plot_units(x_cm,useMicrons)

    for field in fields:
        if ("flash", field) not in ds.field_list and ("flash", field) not in ds.derived_field_list:
            print(f"Skipping {field}: not found")
            continue

        y = squeeze_to_1d(get_field_array(cg, field))

        fig, ax = plt.subplots(figsize=(7, 4))

        ax.plot(x_plot, y, lw=2)

        ax.set_xlabel(xlabel)
        ax.set_ylabel(field)
        ax.set_title(f"{field}, t = {sim_time_ns:.3f} ns")
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        if savePlots:
            out = saveDir / f"{Path(fp).stem}_{field}_1D.png"
            plt.savefig(out, dpi=200)
            print(f"Saved {out}")

        plt.show()


def _auto_downsample_2d(arr, max_pixels=1200):
    n0, n1 = arr.shape
    stride = max(1, int(np.ceil(max(n0, n1) / max_pixels)))
    return arr[::stride, ::stride], stride


def plot_2d_profiles(
    ds,
    fp,
    ftype,
    field,
    useMicrons=True,
    savePlots=False,
    saveDir=Path("."),
    rays=False,
    log2d=True,
):
    """
    Fast plot of one FLASH 2D field + one lineout.

    line_axis = "x" -> fixed x, plot vs y
    line_axis = "y" -> fixed y, plot vs x
    """

    # -----------------------
    # Hardcoded lineout choice
    # -----------------------
    line_axis = "x"       # "x" or "y"
    line_value_um = 0.0

    max_plot_pixels = 1200
    save_dpi = 150

    # -----------------------
    # Load data
    # -----------------------
    cg, dims = get_covering_grid(ds)
    sim_time_ns = get_sim_time_ns(ds)

    arr = get_field_array(cg, ftype, field)
    arr2d = np.asarray(arr[:, :, 0])

    print(f"{field} 2D shape:", arr2d.shape)
    print(f"{field} min/max:", np.nanmin(arr2d), np.nanmax(arr2d))

    # -----------------------
    # Coordinates
    # -----------------------
    x0, x1 = float(ds.domain_left_edge[0]), float(ds.domain_right_edge[0])
    y0, y1 = float(ds.domain_left_edge[1]), float(ds.domain_right_edge[1])

    nx, ny = arr2d.shape

    x_edges_cm = np.linspace(x0, x1, nx + 1)
    y_edges_cm = np.linspace(y0, y1, ny + 1)

    x_cent_cm = 0.5 * (x_edges_cm[:-1] + x_edges_cm[1:])
    y_cent_cm = 0.5 * (y_edges_cm[:-1] + y_edges_cm[1:])

    if useMicrons:
        x_edges = x_edges_cm * 1e4
        y_edges = y_edges_cm * 1e4
        x_cent = x_cent_cm * 1e4
        y_cent = y_cent_cm * 1e4
        line_value = line_value_um
        xlabel = r"$x$ [$\mu$m]"
        ylabel = r"$y$ [$\mu$m]"
        unit_label = r"$\mu$m"
        ray_scale = 1e4
    else:
        x_edges = x_edges_cm
        y_edges = y_edges_cm
        x_cent = x_cent_cm
        y_cent = y_cent_cm
        line_value = line_value_um * 1e-4
        xlabel = r"$x$ [cm]"
        ylabel = r"$y$ [cm]"
        unit_label = "cm"
        ray_scale = 1.0

    # -----------------------
    # Extract lineout first
    # -----------------------
    if line_axis.lower() == "x":
        ix = np.argmin(np.abs(x_cent - line_value))
        line = arr2d[ix, :]
        line_coord = y_cent
        line_location = x_cent[ix]

        line_marker = ("vline", line_location)
        line_xlabel = ylabel
        line_title = rf"Lineout vs $y$ at $x = {line_location:.3f}$ {unit_label}"

    elif line_axis.lower() == "y":
        iy = np.argmin(np.abs(y_cent - line_value))
        line = arr2d[:, iy]
        line_coord = x_cent
        line_location = y_cent[iy]

        line_marker = ("hline", line_location)
        line_xlabel = xlabel
        line_title = rf"Lineout vs $x$ at $y = {line_location:.3f}$ {unit_label}"

    else:
        raise ValueError("line_axis must be either 'x' or 'y'")

    print(f"Lineout axis fixed: {line_axis}")
    print(f"Requested lineout value: {line_value_um:.3f} microns")
    print(f"Nearest lineout location used: {line_location:.4f} {unit_label}")

    # -----------------------
    # 2D display array
    # -----------------------
    arr_img, stride = _auto_downsample_2d(arr2d, max_pixels=max_plot_pixels)

    if stride > 1:
        print(f"Display downsampled by stride {stride}; lineout uses full resolution.")

    finite = arr_img[np.isfinite(arr_img)]
    positive = finite[finite > 0]

    if finite.size == 0:
        raise ValueError(f"{field} has no finite values to plot.")

    if log2d and positive.size > 0:
        vmin = np.nanmin(positive)
        vmax = np.nanmax(finite)
        if vmax <= vmin:
            norm = Normalize(vmin=np.nanmin(finite), vmax=np.nanmax(finite))
        else:
            norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = Normalize(vmin=np.nanmin(finite), vmax=np.nanmax(finite))

    # -----------------------
    # Rays
    # -----------------------
    ray_data = read_flash_rays(fp) if rays else None

    # -----------------------
    # Make figure
    # -----------------------
    fig, (ax2d, axline) = plt.subplots(
        2, 1,
        figsize=(7, 8),
        gridspec_kw={"height_ratios": [3, 1]},
        constrained_layout=True,
    )

    im = ax2d.imshow(
        arr_img.T,
        extent=[x_edges[0], x_edges[-1], y_edges[0], y_edges[-1]],
        origin="lower",
        aspect="equal",
        interpolation="nearest",
        cmap="plasma",
        norm=norm,
    )

    if rays and ray_data is not None:
        plot_rays(
            ax2d,
            ray_data,
            xscale=ray_scale,
            yscale=ray_scale,
            color="w",
            lw=0.8,
        )

    if line_marker[0] == "vline":
        ax2d.axvline(line_marker[1], color="w", ls="--", lw=1.2)
    else:
        ax2d.axhline(line_marker[1], color="w", ls="--", lw=1.2)

    ax2d.set_xlabel(xlabel)
    ax2d.set_ylabel(ylabel)
    ax2d.set_title(f"{field}, t = {sim_time_ns:.3f} ns")

    cbar = fig.colorbar(im, ax=ax2d)
    cbar.set_label(field)

    axline.plot(line_coord, line, lw=2)
    axline.set_xlabel(line_xlabel)
    # -------------- Set x-limits ---------------- #
    axline.set_xlim([200,600])
    # ------------------------------------------- #
    if field == "dens":
        axline.set_ylabel(r"$\rho$ [g/cm$^3$]")
    elif field == "pres":
        axline.set_ylabel(r"$P$ [dyn/cm$^2$]")
    elif field in ["tele", "tion", "trad"]:
        axline.set_ylabel(r"$T$ [eV]")
    elif field in ["velx", "vely"]:
        axline.set_ylabel(r"$v$ [cm/s]")
    elif field == "depo":
        axline.set_ylabel(r"$E_{\rm dep}$ [erg/g]")
    else:
        axline.set_ylabel(field)

    axline.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    axline.yaxis.get_offset_text().set_fontsize(11)

    axline.set_title(line_title)
    axline.grid(True, alpha=0.3)
    axline.set_ylim(0, 1e5)
    # finite_line = line[np.isfinite(line)]
    # if finite_line.size > 0:
    #     ymin = np.nanmin(finite_line)
    #     ymax = np.nanmax(finite_line)
    #     if ymin != ymax:
    #         pad = 0.05 * (ymax - ymin)
    #         axline.set_ylim(ymin - pad, ymax + pad)

    if savePlots:
        saveDir = Path(saveDir)
        saveDir.mkdir(parents=True, exist_ok=True)

        out = saveDir / (
            f"{Path(fp).stem}_{field}_2D_lineout_"
            f"{line_axis}{line_value_um:.1f}um.png"
        )
        plt.savefig(out, dpi=save_dpi, bbox_inches="tight")
        print(f"Saved {out}")

    #plt.show()

    return fig, (ax2d, axline), line_coord, line


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


def plot_rays(ax, ray_data, xscale=1.0, yscale=1.0, color="w", lw=1.0):
    """
    Overlay rays on an existing matplotlib axis.

    xscale, yscale let you convert units.
    Example: FLASH cm -> micron means xscale = 1e4
    """
    tags = ray_data[:, 0]
    unique_tags = np.unique(tags)

    for tag in unique_tags:
        this_ray = ray_data[tags == tag]

        x = this_ray[:, 1] * xscale
        y = this_ray[:, 2] * yscale

        ax.plot(x, y, color=color, lw=lw, alpha=0.8)


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