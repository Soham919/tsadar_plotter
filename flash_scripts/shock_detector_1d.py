
import numpy as np

def detect_simple_density_jump(
    coord, sim_time, dens, tele, tion, velx,
    jump_threshold=5.0e-05,
    post_offset=0.0,
    min_coord=None,
    max_coord=None,
):
    """
    Detect the first density jump scanning right to left.

    All input quantities are in FLASH CGS units,
    except sim_time, which is in ns.

    Parameters
    ----------
    coord : array
        Radial coordinates (cm).
    sim_time : float
        Simulation time (ns).
    dens : array
        Density (g/cm^3).
    tele, tion : arrays
        Electron and ion temperatures (K).
    velx : array
        Radial velocity (cm/s).
    jump_threshold : float
        Minimum density increase (g/cm^3).
    post_offset : float
        Distance behind shock to sample (cm).
    min_coord, max_coord : float or None
        Shock detection bounds (cm).

    Returns
    -------
    result : dict
        Time: ns
        shock_coord: um
        postshock_coord: um
        postshock_vely: cm/s
        postshock_tele: K
        postshock_tion: K
        postshock_dens: g/cm^3
        jump_value: g/cm^3
    """

    if post_offset < 0:
        raise ValueError("post_offset must be >= 0.")

    arrays = [
        np.asarray(a, dtype=float)
        for a in (coord, dens, tele, tion, velx)
    ]

    if any(a.ndim != 1 for a in arrays):
        raise ValueError("All input arrays must be 1D.")

    if len(set(len(a) for a in arrays)) != 1:
        raise ValueError("All input arrays must have equal length.")

    coord, dens, tele, tion, velx = arrays

    # Conversion from cm to micrometers
    CM_TO_UM = 1e4

    result = {
        "found": False,
        "shock_coord": np.nan,
        "Time": sim_time,
        "shock_index": None,
        "jump_value": np.nan,
        "postshock_dens": np.nan,
        "postshock_tele": np.nan,
        "postshock_tion": np.nan,
        "postshock_vely": np.nan,
        "postshock_coord": np.nan,
    }

    # Filter invalid values
    good = np.logical_and.reduce([
        np.isfinite(a) for a in arrays
    ])

    coord, dens, tele, tion, velx = (
        a[good] for a in arrays
    )

    if len(coord) < 2:
        return result

    # Sort from inner to outer radius
    order = np.argsort(coord)

    coord, dens, tele, tion, velx = (
        a[order]
        for a in (coord, dens, tele, tion, velx)
    )

    # Scan from right to left
    for i in range(len(coord) - 2, -1, -1):

        # Shock position in cm (internal calculation)
        shock_coord = 0.5 * (coord[i] + coord[i + 1])

        if min_coord is not None and shock_coord < min_coord:
            continue

        if max_coord is not None and shock_coord > max_coord:
            continue

        jump = dens[i] - dens[i + 1]

        if jump >= jump_threshold:

            # Store shock position in micrometers
            result.update({
                "found": True,
                "shock_coord": float(shock_coord * CM_TO_UM),
                "shock_index": int(i),
                "jump_value": float(jump),
            })

            # Sampling location in cm
            sample_coord = shock_coord - post_offset

            post_indices = np.arange(i + 1)

            if sample_coord < coord[0]:
                return result

            sample_idx = post_indices[
                np.argmin(
                    np.abs(coord[post_indices] - sample_coord)
                )
            ]

            result.update({
                "postshock_dens": float(dens[sample_idx]),
                "postshock_tele": float(tele[sample_idx]),
                "postshock_tion": float(tion[sample_idx]),
                "postshock_vely": float(velx[sample_idx]),

                # Store sampling position in micrometers
                "postshock_coord": float(
                    coord[sample_idx] * CM_TO_UM
                ),
            })

            return result

    return result
