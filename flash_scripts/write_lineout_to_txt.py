import numpy as np
import h5py
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import constants
import yt
import sys
# ----- Mac -----#
#sys.path.append("/Users/soham/Documents/Kinshock/Scripts")
# ----- Windows -----#
sys.path.append(r"\\profiles\Users$\sban\Documents\Scripts")
from flash_scripts.flash_helpers import *


FLASH_UNITS = {
    "r_cm": "cm",
    "dens": "g/cm^3",
    "tele": "K",
    "tion": "K",
    "trad": "K",
    "pres": "erg/cm^3",
    "velx": "cm/s",
    "vely": "cm/s",
    "velr": "cm/s",
    "ye": "dimensionless",
    "cham": "fraction",
    "gas": "fraction",
    "targ": "fraction",
}

#--------- Constants ---------- #
mp = constants.proton_mass # mass of H(kg)
eps = constants.epsilon_0 # epislon naught SI
kb = constants.Boltzmann # Boltzmann constant SI
e = constants.e # charge of electron(Coulomb)
me = constants.electron_mass # mass electron(kg)
mp = constants.proton_mass # mass of H(kg)
pi = constants.pi # pi
c = constants.c

def write_1d_radial_profile(ds, fp, ftype, saveDir=Path("."), fields=None,):
    """
    Extract a full-resolution 1D radial profile from a FLASH dataset
    and write physical quantities to a text file.

    Parameters
    ----------
    ds
        Loaded yt FLASH dataset.

    fp
        Path to the FLASH plotfile.

    ftype
        FLASH field type passed to get_field_array().

    saveDir
        Directory where the text file is written.

    fields : list[str], optional
        FLASH physical fields to export.

        Default:
            dens
            tele
            tion
            trad
            pres
            velx
            ye

    Returns
    -------
    output_path : pathlib.Path
        Path to the generated text file.

    data : dict
        Dictionary containing all exported radial profiles.
    """

    if fields is None:
        fields = [
            "dens",
            "tele",
            "tion",
            "trad",
            "pres",
            "velx",
            "ye"
        ]

    # ----------------------------------------------------------
    # Covering grid
    # ----------------------------------------------------------
    cg, dims = get_covering_grid(ds)

    # Use density to determine the 1D grid size
    dens_arr = np.asarray(
        get_field_array(cg, ftype, "dens")
    ).squeeze()

    if dens_arr.ndim != 1:
        raise ValueError(
            f"Expected a 1D density array, got shape {dens_arr.shape}."
        )

    nr = dens_arr.size

    # ----------------------------------------------------------
    # Cell-centered radial coordinate in cm
    # ----------------------------------------------------------
    r0_cm = float(ds.domain_left_edge[0])
    r1_cm = float(ds.domain_right_edge[0])

    r_edges_cm = np.linspace(
        r0_cm,
        r1_cm,
        nr + 1,
    )

    r_cent_cm = 0.5 * (
        r_edges_cm[:-1] + r_edges_cm[1:]
    )

    # ----------------------------------------------------------
    # Helper for extracting one 1D field
    # ----------------------------------------------------------
    def extract_profile(field_name):

        arr = np.asarray(
            get_field_array(
                cg,
                ftype,
                field_name,
            )
        ).squeeze()

        if arr.ndim != 1:
            raise ValueError(
                f"{field_name} is not 1D after squeeze; "
                f"got shape {arr.shape}."
            )

        if arr.size != nr:
            raise ValueError(
                f"{field_name} has {arr.size} points, "
                f"but density has {nr}."
            )

        return arr

    # ----------------------------------------------------------
    # Extract physical fields
    # ----------------------------------------------------------
    data = {
        "r_cm": r_cent_cm
    }

    for field in fields:
        data[field] = extract_profile(field)

    # In 1D spherical FLASH, velx is radial velocity
    if "velx" in data:
        data["velr"] = data["velx"]
    else:
        data["velr"] = extract_profile("velx")

    # ----------------------------------------------------------
    # Sort by radius
    # ----------------------------------------------------------
    sort_idx = np.argsort(data["r_cm"])

    for key in data:
        data[key] = np.asarray(data[key])[sort_idx]

    # ----------------------------------------------------------
    # Construct output table
    # ----------------------------------------------------------
    column_names = [
        "r_cm",
        *fields,
        "velr",
    ]

    table = np.column_stack(
        [data[name] for name in column_names]
    )

    # ----------------------------------------------------------
    # Save
    # ----------------------------------------------------------
    saveDir = Path(saveDir)

    saveDir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_path = saveDir / (
        f"{Path(fp).stem}_profile_1D.txt"
    )

    header = " ".join(f"{name}[{FLASH_UNITS.get(name, 'unknown')}]" for name in column_names)

    np.savetxt(
        output_path,
        table,
        header=header,
        fmt="%.10e",
    )

    print(f"Saved profile to {output_path}")
    print(f"Columns: {header}")
    print(f"Number of radial points: {nr}")
    print(
        f"Radial range: "
        f"{data['r_cm'][0]:.6e} to "
        f"{data['r_cm'][-1]:.6e} cm"
    )

    return output_path, data


def write_1d_radial_profile_with_species(
    ds,
    fp,
    ftype,
    saveDir=Path("."),
    fields=None,
):
    """
    Extract a full-resolution 1D radial profile from a FLASH dataset
    and write physical quantities + material fractions to a text file.

    Parameters
    ----------
    ds
        Loaded yt FLASH dataset.

    fp
        Path to the FLASH plotfile.

    ftype
        FLASH field type passed to get_field_array().

    saveDir
        Directory where the text file is written.

    fields : list[str], optional
        FLASH physical fields to export.

        Default:
            dens
            tele
            tion
            pres
            velx

    Returns
    -------
    output_path : pathlib.Path
        Path to the generated text file.

    data : dict
        Dictionary containing all exported radial profiles.
    """

    if fields is None:
        fields = [
            "dens",
            "tele",
            "tion",
            "pres",
            "velx",
        ]

    # ----------------------------------------------------------
    # Covering grid
    # ----------------------------------------------------------
    cg, dims = get_covering_grid(ds)

    # Use density to determine the 1D grid size
    dens_arr = np.asarray(
        get_field_array(cg, ftype, "dens")
    ).squeeze()

    if dens_arr.ndim != 1:
        raise ValueError(
            f"Expected a 1D density array, got shape {dens_arr.shape}."
        )

    nr = dens_arr.size

    # ----------------------------------------------------------
    # Cell-centered radial coordinate in cm
    # ----------------------------------------------------------
    r0_cm = float(ds.domain_left_edge[0])
    r1_cm = float(ds.domain_right_edge[0])

    r_edges_cm = np.linspace(
        r0_cm,
        r1_cm,
        nr + 1,
    )

    r_cent_cm = 0.5 * (
        r_edges_cm[:-1] + r_edges_cm[1:]
    )

    # ----------------------------------------------------------
    # Helper for extracting one 1D field
    # ----------------------------------------------------------
    def extract_profile(field_name):

        arr = np.asarray(
            get_field_array(
                cg,
                ftype,
                field_name,
            )
        ).squeeze()

        if arr.ndim != 1:
            raise ValueError(
                f"{field_name} is not 1D after squeeze; "
                f"got shape {arr.shape}."
            )

        if arr.size != nr:
            raise ValueError(
                f"{field_name} has {arr.size} points, "
                f"but density has {nr}."
            )

        return arr

    # ----------------------------------------------------------
    # Extract requested physical fields
    # ----------------------------------------------------------
    data = {
        "r_cm": r_cent_cm
    }

    for field in fields:
        data[field] = extract_profile(field)

    # ----------------------------------------------------------
    # Material fractions
    # ----------------------------------------------------------
    for species in [
        "cham",
        "gas",
        "targ",
    ]:
        data[species] = extract_profile(species)

    # ----------------------------------------------------------
    # Radial velocity
    #
    # In 1D spherical FLASH, velx is the radial velocity.
    # ----------------------------------------------------------
    if "velx" in data:
        data["velr"] = data["velx"]
    else:
        data["velr"] = extract_profile("velx")

    # ----------------------------------------------------------
    # Sort by radius
    # ----------------------------------------------------------
    sort_idx = np.argsort(data["r_cm"])

    for key in data:
        data[key] = np.asarray(data[key])[sort_idx]

    # ----------------------------------------------------------
    # Construct output table
    # ----------------------------------------------------------
    column_names = [
        "r_cm",
        *fields,
        "velr",
        "cham",
        "gas",
        "targ",
    ]

    table = np.column_stack(
        [data[name] for name in column_names]
    )

    # ----------------------------------------------------------
    # Save
    # ----------------------------------------------------------
    saveDir = Path(saveDir)

    saveDir.mkdir(
        parents=True,
        exist_ok=True,
    )

    output_path = saveDir / (
        f"{Path(fp).stem}_profile_1D.txt"
    )

    header = " ".join(f"{name}[{FLASH_UNITS.get(name, 'unknown')}]" for name in column_names)

    np.savetxt(
        output_path,
        table,
        header=header,
        fmt="%.10e",
    )

    print(f"Saved profile to {output_path}")
    print(f"Columns: {header}")
    print(f"Number of radial points: {nr}")
    print(
        f"Radial range: "
        f"{data['r_cm'][0]:.6e} to "
        f"{data['r_cm'][-1]:.6e} cm"
    )

    return output_path, data


def write_2d_lineout_profile(ds, fp, ftype, lineout_dir="x", lineout_pos=0.0, useMicrons=True, saveDir=Path("."), fields=None,):
    """
    Extract a full-resolution 1D lineout from a 2D FLASH dataset
    and write physical quantities + material fractions to a text file.

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
            Hold x fixed at lineout_pos and extract versus y.

        lineout_dir="y":
            Hold y fixed at lineout_pos and extract versus x.

    lineout_pos : float
        Requested fixed-coordinate position.

        In microns when useMicrons=True.
        In cm when useMicrons=False.

    useMicrons : bool
        Controls the units used for lineout_pos only.

        The output coordinate written to file is ALWAYS in cm.

    saveDir
        Directory where the text file is written.

    fields : list[str], optional
        FLASH physical fields to export.

        Default:
            dens
            tele
            tion
            pres
            velx
            vely

    Returns
    -------
    output_path : pathlib.Path
        Path to the generated text file.

    data : dict
        Dictionary containing all exported lineouts.
    """

    lineout_dir = lineout_dir.lower()

    if lineout_dir not in {"x", "y"}:
        raise ValueError(
            "lineout_dir must be either 'x' or 'y'."
        )

    if fields is None:
        fields = [
            "dens",
            "tele",
            "tion",
            "trad",
            "pres",
            "vely",
            "ye"
        ]

    # ----------------------------------------------------------
    # Covering grid
    # ----------------------------------------------------------
    cg, dims = get_covering_grid(ds)

    # Use density only to determine the grid shape
    dens2d = np.asarray(
        get_field_array(cg, ftype, "dens")[:, :, 0]
    )

    nx, ny = dens2d.shape

    # ----------------------------------------------------------
    # Cell-centered coordinates in cm
    # ----------------------------------------------------------
    x0_cm = float(ds.domain_left_edge[0])
    x1_cm = float(ds.domain_right_edge[0])

    y0_cm = float(ds.domain_left_edge[1])
    y1_cm = float(ds.domain_right_edge[1])

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

    x_cent_cm = 0.5 * (
        x_edges_cm[:-1] + x_edges_cm[1:]
    )

    y_cent_cm = 0.5 * (
        y_edges_cm[:-1] + y_edges_cm[1:]
    )

    # ----------------------------------------------------------
    # Convert requested lineout position to cm
    # ----------------------------------------------------------
    if useMicrons:
        requested_pos_cm = lineout_pos * 1.0e-4
    else:
        requested_pos_cm = lineout_pos

    # ----------------------------------------------------------
    # Determine lineout index and coordinate
    # ----------------------------------------------------------
    if lineout_dir == "x":

        # Hold x fixed, lineout varies along y
        line_index = int(
            np.argmin(
                np.abs(x_cent_cm - requested_pos_cm)
            )
        )

        radial_coordinate_cm = y_cent_cm
        actual_position_cm = x_cent_cm[line_index]

    else:

        # Hold y fixed, lineout varies along x
        line_index = int(
            np.argmin(
                np.abs(y_cent_cm - requested_pos_cm)
            )
        )

        radial_coordinate_cm = x_cent_cm
        actual_position_cm = y_cent_cm[line_index]

    print(
        f"Requested fixed-{lineout_dir} position: "
        f"{requested_pos_cm:.6e} cm"
    )

    print(
        f"Nearest cell center used: "
        f"{actual_position_cm:.6e} cm"
    )

    # ----------------------------------------------------------
    # Helper for extracting one field
    # ----------------------------------------------------------
    def extract_line(field_name):

        arr2d = np.asarray(
            get_field_array(
                cg,
                ftype,
                field_name,
            )[:, :, 0]
        )

        if arr2d.shape != (nx, ny):
            raise ValueError(
                f"{field_name} has shape {arr2d.shape}, "
                f"expected {(nx, ny)}."
            )

        if lineout_dir == "x":
            return np.asarray(
                arr2d[line_index, :]
            )

        else:
            return np.asarray(
                arr2d[:, line_index]
            )

    # ----------------------------------------------------------
    # Extract requested physical fields
    # ----------------------------------------------------------
    data = {
        "r_cm": np.asarray(radial_coordinate_cm)
    }

    for field in fields:
        data[field] = extract_line(field)

    # ----------------------------------------------------------
    # Material fractions
    # ----------------------------------------------------------
    for species in [
        "cham",
        "gas",
        "targ",
    ]:
        data[species] = extract_line(species)

    data["velr"] = data["vely"]

    # ----------------------------------------------------------
    # Sort by coordinate
    # ----------------------------------------------------------
    sort_idx = np.argsort(data["r_cm"])

    for key in data:
        data[key] = np.asarray(data[key])[sort_idx]

    # ----------------------------------------------------------
    # Construct output table
    # ----------------------------------------------------------
    column_names = [
        "r_cm",
        *fields,
        "velr",
        "cham",
        "gas",
        "targ",
    ]

    table = np.column_stack(
        [data[name] for name in column_names]
    )

    # ----------------------------------------------------------
    # Save
    # ----------------------------------------------------------
    saveDir = Path(saveDir)

    saveDir.mkdir(
        parents=True,
        exist_ok=True,
    )

    unit_suffix = (
        "um" if useMicrons else "cm"
    )

    output_path = saveDir / (
        f"{Path(fp).stem}_profile_"
        f"{lineout_dir}_{lineout_pos:g}{unit_suffix}.txt"
    )

    header = " ".join(f"{name}[{FLASH_UNITS.get(name, 'unknown')}]" for name in column_names)

    np.savetxt(
        output_path,
        table,
        header=header,
        fmt="%.10e",
    )

    print(f"Saved profile to {output_path}")
    print(f"Columns: {header}")
    print(f"Number of radial points: {len(radial_coordinate_cm)}")

    return output_path, data

if __name__ == "__main__":
    runDir = Path(r"C:\Simulation_data\FLASH\1D_spherical\iFP_bubble_T_7.5e10")  # directory containing the FLASH plot files

    # Axis labels
    useMicrons = True

    # Save plots?
    savePlots = True
    saveDir = runDir / "plots"
    saveDir.mkdir(parents=True, exist_ok=True)

    target_time_ns = 1  # ns
    files_info = find_nearest_flash_file_from_table(target_time_ns, runDir, time_table="flash_plotfile_times.txt")
               
    filename, file_num, matched_time_s, matched_time_ns = files_info

    fp = filename

    print("\n" + "=" * 70)
    print(f"Loading file number {file_num}")
    print(f"File: {fp}")
    print(f"Matched time = {matched_time_ns:.4f} ns")
    print("=" * 70)

    ds = yt.load(fp)

    #output_path, data = write_2d_lineout_profile(ds, fp, ftype = "flash", lineout_dir="x", lineout_pos= 0, useMicrons=useMicrons, saveDir=saveDir)
    output_path, data = write_1d_radial_profile(ds, fp, ftype = "flash", saveDir=saveDir)

    ## Optional check for profiles you want import 
    plt.plot(data['r_cm'], data['tion'])
    plt.plot(data['r_cm'], data['tele'])
    plt.plot(data['r_cm'], data['trad'])
    #plt.plot()

    plt.show()

