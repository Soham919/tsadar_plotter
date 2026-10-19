import numpy as np
from scipy.interpolate import interp1d


def intext_1D(x, y, kind="linear"):
    """
    Create a function y(x) from tabulated 1D data.

    Parameters
    ----------
    x : array_like
        Argument values, e.g. temperature [eV].

    y : array_like
        Function values, e.g. specific internal energy [erg/g].

    kind : str, optional
        Interpolation type. Default is "linear".

    Returns
    -------
    f : callable
        Function that interpolates inside the tabulated range
        and extrapolates outside it.
    """

    x = np.asarray(x)
    y = np.asarray(y)

    if x.ndim != 1 or y.ndim != 1:
        raise ValueError("x and y must both be 1D arrays.")

    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")

    # Sort in increasing x, just in case
    idx = np.argsort(x)
    x = x[idx]
    y = y[idx]

    f = interp1d(
        x,
        y,
        kind=kind,
        bounds_error=False,
        fill_value="extrapolate"
    )

    return f