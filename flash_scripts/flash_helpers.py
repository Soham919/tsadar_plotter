import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, Normalize
from pathlib import Path
import yt

# ============================================================
# Helper functions
# ============================================================

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
    Returns max-level covering grid and dimensions.
    """
    level = ds.index.max_level
    dims = ds.domain_dimensions * (2 ** level)

    print("max level:", level)
    print("domain_dimensions:", ds.domain_dimensions)
    print("covering grid dims:", dims)

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
