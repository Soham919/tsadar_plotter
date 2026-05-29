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


def plot_2d_field(ds, ftype, field, useMicrons, savePlots, saveDir, fp, rays = False):
    """
    Plot 2D pcolormesh from FLASH output.
    """
    cg, dims = get_covering_grid(ds)
    sim_time_ns = get_sim_time_ns(ds)
    arr = get_field_array(cg, ftype, field)

    print(f"{field} raw shape from covering_grid:", arr.shape)

    # Pick z-slice explicitly
    arr2d = arr[:, :, 0]

    print(f"{field} 2D shape:", arr2d.shape)
    print(f"{field} min/max:", np.nanmin(arr2d), np.nanmax(arr2d))

    # Domain edges
    x0, x1 = float(ds.domain_left_edge[0]), float(ds.domain_right_edge[0])
    y0, y1 = float(ds.domain_left_edge[1]), float(ds.domain_right_edge[1])

    nx, ny = arr2d.shape

    # Cell edges for pcolormesh
    x_edges_cm = np.linspace(x0, x1, nx + 1)
    y_edges_cm = np.linspace(y0, y1, ny + 1)

    if useMicrons:
        x_edges = x_edges_cm * 1e4
        y_edges = y_edges_cm * 1e4
        xlabel = r"$x$ [$\mu$m]"
        ylabel = r"$y$ [$\mu$m]"
    else:
        x_edges = x_edges_cm
        y_edges = y_edges_cm
        xlabel = r"$x$ [cm]"
        ylabel = r"$y$ [cm]"

    X, Y = np.meshgrid(x_edges, y_edges, indexing="ij")

    positive = arr2d[arr2d > 0]

    if rays == True:
        ray_data = read_flash_rays(fp)
    
    fig, ax = plt.subplots(figsize=(6, 6))

    if positive.size > 0:
        norm = LogNorm(vmin=np.nanmin(positive), vmax=np.nanmax(arr2d))
    else:
        print(f"Warning: {field} has no positive values. Using linear scale.")
        norm = Normalize(vmin=np.nanmin(arr2d), vmax=np.nanmax(arr2d))

    pcm = ax.pcolormesh(
        X,
        Y,
        arr2d,
        shading="auto",
        cmap="plasma",
        norm=norm,
    )

    if rays == True :
        # overlay rays
        plot_rays(ax, ray_data, xscale=1e4, yscale=1e4, color="w", lw=0.8)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{field}, t = {sim_time_ns:.3f} ns")
    ax.set_aspect("equal")

    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label(field)

    plt.tight_layout()

    if savePlots:
        out = saveDir / f"{Path(fp).stem}_{field}_2D.png"
        plt.savefig(out, dpi=200)
        print(f"Saved {out}")

    plt.show()


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

