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
    """
    Downsample a 2D array for display only.
    The full-resolution array remains available for lineouts.
    """
    n0, n1 = arr.shape
    stride = max(1, int(np.ceil(max(n0, n1) / max_pixels)),)
    return arr[::stride, ::stride], stride


def plot_2d_with_lineout(
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
        cmap="plasma",
        norm = norm,
        rasterized=True,
    )

    ax_2d.set_xlabel(image_xlabel)
    ax_2d.set_ylabel(image_ylabel)
    ax_2d.set_title(f"{field}, t = {sim_time_ns:.3f} ns")

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

    # Match the horizontal ranges of the two lineout panels.
    line_min = line_coordinate[0]
    line_max = line_coordinate[-1]
    ax_line.set_xlim(line_min,line_max,)
    ax_fraction.set_xlim(line_min,line_max,)

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

    return (fig,(ax_line, ax_2d, ax_fraction),line_data,)


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